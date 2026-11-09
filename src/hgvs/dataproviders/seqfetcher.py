"""provides sequencing fetching from NCBI and Ensembl"""

import logging
import os
import threading
from collections import OrderedDict

import bioutils.seqfetcher

import hgvs

from ..exceptions import HGVSDataNotAvailableError

_logger = logging.getLogger(__name__)


class _Window:
    """A stretch of one sequence held in memory: seq is sequence[start:end]."""

    __slots__ = ("end", "hits", "known_length", "seq", "size", "start")

    def __init__(self, start, seq, size, requested_end):
        self.start = start
        self.end = start + len(seq)
        self.seq = seq
        self.size = size  # read-ahead size in force when this window was fetched
        self.hits = 0  # requests served from this window
        # A fetch that came back short ran into the end of the sequence, so the
        # sequence's length is now known and requests past it need not be windowed.
        self.known_length = self.end if self.end < requested_end else None

    def covers(self, start_i, end_i):
        return self.start <= start_i and end_i <= self.end


class ReadAheadCache:
    """Serves small sequence requests from larger windows fetched from a backend.

    hgvs asks for sequence in tiny pieces: one base to check a variant's reference,
    about 20 bases for the normalizer's shuffling window. Against a remote sequence
    source (SeqRepo REST, NCBI eutils) each piece is a round trip, and those round
    trips dominate the time to map a variant. When requests walk along a sequence,
    as they do for variants sorted by position (a VCF, for example), one fetch of a
    window of several kilobases serves the following hundreds of requests from memory.

    One window is kept per accession, for the most recently used accessions. The
    window size adapts to the access pattern: it doubles after a window that served at
    least one request from memory and halves after one that served none, between
    ``min_size`` and ``max_size``. Sequential access soon reaches ``max_size`` while
    random access stays at ``min_size``, which costs about the same as fetching the
    bases actually asked for. A request walking backwards (descending positions, as
    for a minus-strand transcript's variants in c. order) gets its window placed
    before it.

    The backend's own behaviour is kept for anything unusual: a request that reaches
    past the end of the sequence, a request larger than ``max_size``, or a window
    fetch that fails is passed to the backend as-is, so errors and truncation are
    exactly what the backend would have produced.

    :param fetcher: callable ``(ac, start_i, end_i) -> str`` returning
        ``sequence[start_i:end_i]``, with ``end_i=None`` meaning the end of the sequence
    """

    def __init__(self, fetcher, min_size=4096, max_size=1048576, max_accessions=64):
        if not 0 < min_size <= max_size:
            msg = f"read-ahead sizes must satisfy 0 < min_size <= max_size ({min_size}, {max_size})"
            raise ValueError(msg)
        self.fetcher = fetcher
        self.min_size = min_size
        self.max_size = max_size
        self.max_accessions = max_accessions
        self._windows = OrderedDict()  # ac -> _Window, least recently used first
        self._lock = threading.Lock()

    @classmethod
    def from_config(cls, fetcher, config=None):
        """Build a cache from the ``[seqfetcher]`` section of ``hgvs.global_config``
        (or ``config``). Returns None when ``readahead`` is turned off there."""
        cfg = config if config is not None else hgvs.global_config.seqfetcher
        if not cfg.readahead:
            return None
        return cls(
            fetcher,
            min_size=cfg.readahead_min_size,
            max_size=cfg.readahead_max_size,
            max_accessions=cfg.readahead_max_accessions,
        )

    def fetch(self, ac, start_i=None, end_i=None):
        """Return ``sequence[start_i:end_i]`` for ``ac``, from a window when one covers it."""
        if start_i is None:
            start_i = 0
        if end_i is None:
            return self._fetch_to_end(ac, start_i)
        if start_i < 0 or end_i < start_i or end_i - start_i > self.max_size:
            return self.fetcher(ac, start_i, end_i)
        with self._lock:
            window = self._windows.get(ac)
            if window is not None and window.covers(start_i, end_i):
                window.hits += 1
                self._windows.move_to_end(ac)
                return window.seq[start_i - window.start : end_i - window.start]
        return self._fetch_window(ac, start_i, end_i, window)

    def _fetch_to_end(self, ac, start_i):
        # Whole sequence (or its tail): fetch directly and keep it when it's small
        # enough. Later requests inside it (common for transcripts) are then free.
        seq = self.fetcher(ac, start_i, None)
        if len(seq) <= self.max_size:
            self._store(ac, _Window(start_i, seq, self.min_size, float("inf")))
        return seq

    def _fetch_window(self, ac, start_i, end_i, previous):
        """Fetch a window around ``[start_i, end_i)``, keep it, and return the slice.

        ``previous`` is the accession's current window (which did not cover the
        request) or None. It sets the size and direction of the new window.
        """
        if (
            previous is not None
            and previous.known_length is not None
            and end_i > previous.known_length
        ):
            # Past the end of the sequence: the backend decides what that means.
            return self.fetcher(ac, start_i, end_i)
        size = self._next_size(previous)
        # The window starts a little before the request (or, walking backwards, ends a
        # little after it): having checked one base, hgvs next asks for the bases
        # around it, and that should be a hit too.
        margin = size // 16
        if previous is not None and start_i < previous.start:
            w_end = end_i + margin
            w_start = min(start_i, max(0, w_end - size))
            w_end = max(w_end, w_start + size)
        else:
            w_start = max(0, start_i - margin)
            w_end = max(w_start + size, end_i)
        try:
            seq = self.fetcher(ac, w_start, w_end)
        except Exception as ex:
            _logger.debug("Read-ahead fetch of %s[%d:%d] failed (%s)", ac, w_start, w_end, ex)
            return self.fetcher(ac, start_i, end_i)
        window = _Window(w_start, seq, size, w_end)
        self._store(ac, window)
        if window.covers(start_i, end_i):
            return seq[start_i - w_start : end_i - w_start]
        # The sequence ended before the request did: ask the backend for exactly this
        # range, so that out-of-bounds requests behave as they would without read-ahead.
        return self.fetcher(ac, start_i, end_i)

    def prefetch(self, ac, start_i, end_i):
        """Fetch ``sequence[start_i:end_i]`` in one request and keep it for later requests.

        For callers that know the region they are about to work in. There is no size
        limit; the caller is trusted to ask for something sensible.
        """
        seq = self.fetcher(ac, start_i, end_i)
        self._store(ac, _Window(start_i, seq, self.min_size, end_i))

    def clear(self):
        with self._lock:
            self._windows.clear()

    def _next_size(self, previous):
        if previous is None:
            return self.min_size
        if previous.hits:
            return min(previous.size * 2, self.max_size)
        return max(previous.size // 2, self.min_size)

    def _store(self, ac, window):
        with self._lock:
            self._windows[ac] = window
            self._windows.move_to_end(ac)
            while len(self._windows) > self.max_accessions:
                self._windows.popitem(last=False)


class ReadAheadSeqFetcher:
    """Put a :class:`ReadAheadCache` in front of any sequence fetcher object.

    For data providers that bring their own fetcher instead of :class:`SeqFetcher`
    (cdot's chained, prefix and FASTA fetchers, for example). The wrapped object needs
    ``fetch_seq(ac, start_i, end_i)``; ``source`` and ``set_data_provider`` are passed
    through when it has them, so the wrapper drops in wherever the original did::

        hdp = JSONDataProvider(..., seqfetcher=ReadAheadSeqFetcher(my_fetcher))

    Window sizes come from ``hgvs.global_config.seqfetcher`` unless given here. With
    ``readahead`` turned off in the config, requests pass straight through.
    """

    def __init__(self, seqfetcher, min_size=None, max_size=None, max_accessions=None):
        self.seqfetcher = seqfetcher
        cfg = hgvs.global_config.seqfetcher
        self._readahead = None
        if cfg.readahead:
            self._readahead = ReadAheadCache(
                seqfetcher.fetch_seq,
                min_size=min_size if min_size is not None else cfg.readahead_min_size,
                max_size=max_size if max_size is not None else cfg.readahead_max_size,
                max_accessions=(
                    max_accessions if max_accessions is not None else cfg.readahead_max_accessions
                ),
            )

    @property
    def source(self):
        return getattr(self.seqfetcher, "source", type(self.seqfetcher).__name__)

    def set_data_provider(self, hdp):
        set_data_provider = getattr(self.seqfetcher, "set_data_provider", None)
        if set_data_provider is not None:
            set_data_provider(hdp)

    def fetch_seq(self, ac, start_i=None, end_i=None):
        if self._readahead is None:
            return self.seqfetcher.fetch_seq(ac, start_i, end_i)
        return self._readahead.fetch(ac, start_i, end_i)

    def prefetch(self, ac, start_i, end_i):
        """Fetch ``ac[start_i:end_i]`` once so later requests inside it are served
        from memory. Does nothing when read-ahead is disabled."""
        if self._readahead is not None:
            self._readahead.prefetch(ac, start_i, end_i)


class SeqFetcher:
    """This class is intended primarily as a mixin for HGVS data providers
    that doen't otherwise have access to sequence data.  It uses the
    fetch_seq() function in this module to fetch sequences from
    several sources; see that function for details.

    Requests are served through a :class:`ReadAheadCache` unless
    ``hgvs.global_config.seqfetcher.readahead`` is False; see that class for why.

    >> sf = SeqFetcher()

    >> sf.fetch_seq('NP_056374.2',0,10)
    'MESRETLSSS'

    """

    def __init__(self):
        # If HGVS_SEQREPO_DIR is defined, we use seqrepo for *all* sequences.
        # If HGVS_SEQREPO_URL is defined, use instance of seqrepo-rest-service for *all* sequences.
        #   (see https://github.com/biocommons/seqrepo-rest-service for more info)
        # Otherwise, we fall back to remote sequence fetching
        seqrepo_dir = os.environ.get("HGVS_SEQREPO_DIR")
        seqrepo_url = os.environ.get("HGVS_SEQREPO_URL")
        if seqrepo_dir:
            from biocommons.seqrepo import SeqRepo

            self.sr = SeqRepo(seqrepo_dir)

            def _fetch_seq_seqrepo(ac, start_i=None, end_i=None):
                return self.sr.fetch(ac, start_i, end_i)

            self.fetcher = _fetch_seq_seqrepo
            self.source = f"SeqRepo ({seqrepo_dir})"
        elif seqrepo_url:
            from biocommons.seqrepo.dataproxy import SeqRepoRESTDataProxy

            self.sr = SeqRepoRESTDataProxy(seqrepo_url)
            self.fetcher = lambda ac, start_i=None, end_i=None: self.sr.get_sequence(
                ac, start_i, end_i
            )
            self.source = f"SeqRepo REST ({seqrepo_url})"
        else:
            self.sr = None
            self.fetcher = bioutils.seqfetcher.fetch_seq
            self.source = "bioutils.seqfetcher (network fetching)"

        self._readahead = ReadAheadCache.from_config(self.fetcher)
        if self._readahead is not None:
            _logger.info(
                "Fetching sequences with %s, read-ahead %d-%d bp",
                self.source,
                self._readahead.min_size,
                self._readahead.max_size,
            )
        else:
            _logger.info("Fetching sequences with %s", self.source)

    def fetch_seq(self, ac, start_i=None, end_i=None):
        fetch = self._readahead.fetch if self._readahead is not None else self.fetcher
        try:
            return fetch(ac, start_i, end_i)
        except Exception as ex:
            msg = f"Failed to fetch {ac} from {self.source} ({ex})"
            raise HGVSDataNotAvailableError(msg) from ex

    def prefetch(self, ac, start_i, end_i):
        """Fetch ``ac[start_i:end_i]`` in one request so later requests within it are
        served from memory. Useful when the variants to come are known to lie in a
        region but are not sorted. Does nothing when read-ahead is disabled.
        """
        if self._readahead is None:
            return
        try:
            self._readahead.prefetch(ac, start_i, end_i)
        except Exception as ex:
            msg = f"Failed to fetch {ac} from {self.source} ({ex})"
            raise HGVSDataNotAvailableError(msg) from ex


# <LICENSE>
# Copyright 2018 HGVS Contributors (https://github.com/biocommons/hgvs)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
