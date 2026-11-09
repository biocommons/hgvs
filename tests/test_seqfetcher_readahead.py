"""ReadAheadCache: serve hgvs' tiny sequence requests from larger windows."""

import random

import pytest

import hgvs
from hgvs.dataproviders.seqfetcher import ReadAheadCache, ReadAheadSeqFetcher, SeqFetcher
from hgvs.exceptions import HGVSDataNotAvailableError


class FakeBackend:
    """In-memory sequences with slice semantics, recording every request."""

    def __init__(self, **seqs):
        self.seqs = seqs
        self.calls = []

    def __call__(self, ac, start_i=None, end_i=None):
        self.calls.append((ac, start_i, end_i))
        if ac not in self.seqs:
            msg = f"no such sequence {ac}"
            raise KeyError(msg)
        if start_i is not None and start_i > len(self.seqs[ac]):
            msg = f"start {start_i} beyond end of {ac}"
            raise ValueError(msg)
        return self.seqs[ac][start_i:end_i]


def _random_seq(n, seed):
    return "".join(random.Random(seed).choices("ACGT", k=n))  # noqa: S311


CHR = _random_seq(50_000, 1)
TX = _random_seq(918, 2)


@pytest.fixture
def backend():
    return FakeBackend(CHR=CHR, TX=TX)


@pytest.fixture
def cache(backend):
    return ReadAheadCache(backend, min_size=1024, max_size=8192, max_accessions=3)


def test_sorted_requests_share_one_window(backend, cache):
    # 1-base checks and 20-base normalizer windows walking along the sequence
    for pos in range(100, 1000, 7):
        assert cache.fetch("CHR", pos, pos + 1) == CHR[pos : pos + 1]
        assert cache.fetch("CHR", pos - 10, pos + 11) == CHR[pos - 10 : pos + 11]
    assert backend.calls == [("CHR", 100 - 64, 100 - 64 + 1024)]  # starts a little early


def test_window_grows_while_requests_keep_landing_in_it(backend, cache):
    pos = 0
    while pos < 40_000:
        assert cache.fetch("CHR", pos, pos + 1) == CHR[pos : pos + 1]
        pos += 100
    sizes = [e - s for _, s, e in backend.calls]
    assert sizes[:5] == [1024, 2048, 4096, 8192, 8192]  # doubles up to max_size, then stays
    assert all(e <= 50_000 + 8192 for _, _, e in backend.calls)
    # ~40 kb of sequential 1-base requests cost a handful of fetches, not 400
    assert len(backend.calls) <= 10


def test_window_stays_small_for_scattered_access(backend, cache):
    for pos in range(39_000, 0, -3000):  # each request outside the previous window
        assert cache.fetch("CHR", pos, pos + 1) == CHR[pos : pos + 1]
    sizes = {e - s for _, s, e in backend.calls}
    assert sizes == {1024}  # never grew: no window served a second request
    assert len(backend.calls) == 13


def test_window_placed_before_request_when_walking_backwards(backend, cache):
    cache.fetch("CHR", 30_000, 30_001)
    cache.fetch("CHR", 30_050, 30_051)  # a hit, so the next window doubles
    assert cache.fetch("CHR", 20_000, 20_001) == CHR[20_000:20_001]
    assert backend.calls[-1] == ("CHR", 20_001 + 128 - 2048, 20_001 + 128)
    assert cache.fetch("CHR", 19_000, 19_001) == CHR[19_000:19_001]  # served from that window
    assert len(backend.calls) == 2


def test_window_does_not_start_before_zero(backend, cache):
    cache.fetch("CHR", 5000, 5001)
    assert cache.fetch("CHR", 10, 11) == CHR[10:11]
    assert backend.calls[-1] == ("CHR", 0, 1024)


def test_results_match_backend_for_every_request_shape(cache):
    rnd = random.Random(4)  # noqa: S311
    for _ in range(500):
        start = rnd.randrange(0, 50_000)
        end = min(50_000, start + rnd.choice([1, 1, 1, 21, 21, 300, 5000, 20_000]))
        assert cache.fetch("CHR", start, end) == CHR[start:end]


def test_request_past_end_of_sequence_is_answered_by_backend(backend, cache):
    assert cache.fetch("TX", 900, 918) == TX[900:918]
    assert backend.calls == [("TX", 836, 1860)]  # window ran into the end of TX
    backend.calls.clear()
    # Out of bounds: exactly what the backend returns, and no window fetch first
    assert cache.fetch("TX", 910, 930) == TX[910:930]
    assert cache.fetch("TX", 918, 919) == ""
    assert backend.calls == [("TX", 910, 930), ("TX", 918, 919)]


def test_short_fetch_on_first_window_still_answers_request_from_backend(backend, cache):
    # the request itself reaches past the end of the sequence
    assert cache.fetch("TX", 910, 930) == TX[910:930]
    assert backend.calls == [("TX", 846, 1870), ("TX", 910, 930)]


def test_whole_sequence_fetch_is_kept_when_small(backend, cache):
    assert cache.fetch("TX") == TX
    assert cache.fetch("TX", 100, 101) == TX[100:101]
    assert backend.calls == [("TX", 0, None)]


def test_whole_sequence_fetch_is_not_kept_when_large(backend, cache):
    assert cache.fetch("CHR") == CHR
    cache.fetch("CHR", 100, 101)
    assert backend.calls == [("CHR", 0, None), ("CHR", 36, 1060)]


def test_oversized_request_bypasses_the_cache(backend, cache):
    assert cache.fetch("CHR", 0, 20_000) == CHR[:20_000]
    assert backend.calls == [("CHR", 0, 20_000)]
    cache.fetch("CHR", 5, 6)
    assert backend.calls[-1] == ("CHR", 0, 1024)


def test_failed_window_fetch_falls_back_to_exact_request(cache):
    class Flaky(FakeBackend):
        def __call__(self, ac, start_i=None, end_i=None):
            if end_i is not None and end_i - start_i > 100:
                self.calls.append((ac, start_i, end_i))
                msg = "range too large"
                raise RuntimeError(msg)
            return super().__call__(ac, start_i, end_i)

    backend = Flaky(CHR=CHR)
    cache.fetcher = backend
    assert cache.fetch("CHR", 10, 31) == CHR[10:31]
    assert backend.calls == [("CHR", 0, 1024), ("CHR", 10, 31)]


def test_backend_errors_for_the_exact_request_propagate(cache):
    with pytest.raises(KeyError):
        cache.fetch("NOPE", 0, 1)


def test_least_recently_used_accessions_are_dropped(cache):
    backend = FakeBackend(A=CHR, B=CHR, C=CHR, D=CHR)
    cache.fetcher = backend
    for ac in "ABCD":
        cache.fetch(ac, 0, 1)
    cache.fetch("B", 1, 2)  # hit: B stays warm
    assert list(cache._windows) == ["C", "D", "B"]
    cache.fetch("A", 1, 2)  # A was evicted: fetched again
    assert backend.calls[-1] == ("A", 0, 1024)


def test_prefetch_warms_a_region(backend, cache):
    cache.prefetch("CHR", 1000, 30_000)  # bigger than max_size is fine when explicit
    for pos in range(1000, 30_000, 1234):
        assert cache.fetch("CHR", pos, pos + 1) == CHR[pos : pos + 1]
    assert backend.calls == [("CHR", 1000, 30_000)]


def test_invalid_sizes_are_rejected(backend):
    with pytest.raises(ValueError, match="min_size <= max_size"):
        ReadAheadCache(backend, min_size=10, max_size=5)


@pytest.fixture
def network_free_seqfetcher(monkeypatch, backend):
    monkeypatch.setenv("HGVS_SEQREPO_DIR", "")
    monkeypatch.setenv("HGVS_SEQREPO_URL", "")
    monkeypatch.setattr("bioutils.seqfetcher.fetch_seq", backend)
    return backend


def test_seqfetcher_uses_readahead_by_default(network_free_seqfetcher, monkeypatch):
    monkeypatch.setattr(hgvs.global_config.seqfetcher, "readahead_min_size", 1024)
    sf = SeqFetcher()
    assert sf.fetcher is network_free_seqfetcher  # the raw backend is still exposed
    assert sf.fetch_seq("CHR", 10, 11) == CHR[10:11]
    assert sf.fetch_seq("CHR", 500, 521) == CHR[500:521]
    assert network_free_seqfetcher.calls == [("CHR", 0, 1024)]
    sf.prefetch("CHR", 20_000, 21_000)
    assert sf.fetch_seq("CHR", 20_500, 20_501) == CHR[20_500:20_501]
    assert network_free_seqfetcher.calls[-1] == ("CHR", 20_000, 21_000)
    with pytest.raises(HGVSDataNotAvailableError, match="Failed to fetch NOPE"):
        sf.fetch_seq("NOPE", 0, 1)
    with pytest.raises(HGVSDataNotAvailableError, match="Failed to fetch NOPE"):
        sf.prefetch("NOPE", 0, 1)


def test_seqfetcher_readahead_can_be_disabled(network_free_seqfetcher, monkeypatch):
    monkeypatch.setattr(hgvs.global_config.seqfetcher, "readahead", False)
    sf = SeqFetcher()
    assert sf._readahead is None
    assert sf.fetch_seq("CHR", 10, 11) == CHR[10:11]
    sf.prefetch("CHR", 0, 100)  # no-op
    assert network_free_seqfetcher.calls == [("CHR", 10, 11)]


class FakeCdotStyleFetcher:
    """What cdot passes as ``seqfetcher=``: fetch_seq, source, set_data_provider."""

    source = "fake"

    def __init__(self, backend):
        self.backend = backend
        self.hdp = None

    def set_data_provider(self, hdp):
        self.hdp = hdp

    def fetch_seq(self, ac, start_i=None, end_i=None):
        return self.backend(ac, start_i, end_i)


def test_readahead_seqfetcher_wraps_any_fetch_seq_object(backend, monkeypatch):
    monkeypatch.setattr(hgvs.global_config.seqfetcher, "readahead_min_size", 1024)
    inner = FakeCdotStyleFetcher(backend)
    sf = ReadAheadSeqFetcher(inner)
    assert sf.source == "fake"
    sf.set_data_provider("hdp")
    assert inner.hdp == "hdp"
    assert sf.fetch_seq("CHR", 10, 11) == CHR[10:11]
    assert sf.fetch_seq("CHR", 500, 521) == CHR[500:521]
    assert backend.calls == [("CHR", 0, 1024)]
    sf.prefetch("TX", 0, 918)
    assert sf.fetch_seq("TX", 5, 6) == TX[5:6]
    assert backend.calls[-1] == ("TX", 0, 918)


def test_readahead_seqfetcher_explicit_sizes_and_bare_objects(backend):
    class Bare:
        def fetch_seq(self, ac, start_i=None, end_i=None):
            return backend(ac, start_i, end_i)

    sf = ReadAheadSeqFetcher(Bare(), min_size=256, max_size=512)
    assert sf.source == "Bare"
    sf.set_data_provider(None)  # no-op when the wrapped object has none
    assert sf.fetch_seq("CHR", 100, 101) == CHR[100:101]
    assert backend.calls == [("CHR", 100 - 16, 100 - 16 + 256)]


def test_readahead_seqfetcher_passes_through_when_disabled(backend, monkeypatch):
    monkeypatch.setattr(hgvs.global_config.seqfetcher, "readahead", False)
    sf = ReadAheadSeqFetcher(FakeCdotStyleFetcher(backend))
    assert sf.fetch_seq("CHR", 10, 11) == CHR[10:11]
    sf.prefetch("CHR", 0, 100)
    assert backend.calls == [("CHR", 10, 11)]
