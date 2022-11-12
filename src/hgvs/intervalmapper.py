# -*- coding: utf-8 -*-
"""
Mapping intervals between pairs of congruent segments

The IntervalMapper class is at the heart of mapping between aligned sequences.  An instance
of :class:`uta.tools.intervalmapper.IntervalMapper` is constructed with an ordered list of
:class:`uta.tools.intervalmapper.IntervalPair` instances, each of which consists of two
:class:`uta.tools.intervalmapper.Interval` instances.  The IntervalMapper class is unaware
of strand/orientation; that issue is handled by the
:class:`uta.tools.transcriptmapper.TranscriptMapper` class.

NOTE: Mapping at the boundaries around indels requires a choice.  If seq B
has an insertion relative to seq A, then mapping coordinate at the
boundaries can either be minimal or maximal for both the start and
end. Consider this alignment::

        0         15   20         35         50
        |==========|====|==========|==========|
        |          | __/        __/|          |
        |          |/          /   |          |
        |==========|==========|====|==========|
        0         15         30   35         50
            15M   5D   15M      5I      15M  

  segment 1: [ 0,15] ~ [ 0,15]
  segment 2: [15,20] ~ [15,15]
  segment 3: [20,35] ~ [15,30]
  segment 4: [35,35] ~ [30,35]
  segment 5: [35,50] ~ [35,50]

and these intervals around reference position 35::

  interval 1: 34,36   -> 29,36 (no ambiguity)
  interval 2: 35,35   -> 30,35 (reasonable)
  interval 3: 34,35   -> 29,30 (minimal) or 29,35 (maximal)
  interval 4: 35,36   -> 35,36 (minimal) or 30,36 (maximal)

So, for interval 3, end_i=35 matches segment 3 and segment 4.  Analagously
for interval 4, start_i=35 matches segment 4 and segment 5.

Currently, this code matches an interval <start_i,end_i> using the maximal
start_i and minimal end_i.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import re

from hgvs.exceptions import HGVSInvalidIntervalError
from six.moves import range

_logger = logging.getLogger(__name__)
_logger.warning("This module is deprecated and will be removed in a future release")


# N.B. This Interval is internal to intervalmapper.py. It is *NOT* the
# same as the Interval defined in location.py.
class Interval(object):
    """Represents a segment of a sequence in interbase
    coordinates (0-based, right-open).
    """
    __slots__ = ("start_i", "end_i")

    def __init__(self, start_i, end_i):
        if not (start_i <= end_i):
            raise HGVSInvalidIntervalError("start_i must be less than or equal to end_i")
        self.start_i = start_i
        self.end_i = end_i

    @property
    def len(self):
        return self.end_i - self.start_i

    def __repr__(self):
        return "{self.__class__.__name__}(start_i={self.start_i},end_i={self.end_i})".format(
            self=self)


class IntervalPair(object):
    """Represents a match, insertion, or deletion segment of an
    alignment. If a match, the lengths must be equal; if an insertion or
    deletion, the length of the ref or tgt must be zero respectively."""
    __slots__ = ("ref", "tgt")

    def __init__(self, ref, tgt):
        if not ((ref.len == tgt.len) or (ref.len == 0 and tgt.len != 0) or
                (ref.len != 0 and tgt.len == 0)):
            raise HGVSInvalidIntervalError(
                "IntervalPair doesn't represent a match, insertion, or deletion")
        self.ref = ref
        self.tgt = tgt

    def __repr__(self):
        return "{self.__class__.__name__}(ref={self.ref},tgt={self.tgt})".format(self=self)


class IntervalMapper(object):
    """Provides mapping between sequence coordinates according to an
    ordered set of IntervalPairs."""
    __slots__ = ("interval_pairs", "ref_intervals", "tgt_intervals", "ref_len", "tgt_len")

    def __init__(self, interval_pairs):
        """
        :param interval_pairs: an ordered list of IntervalPair instances
        :type interval_pairs: list (of IntervalPair instances).
        :returns: an IntervalMapper instance
        """

        def _validate_intervals(ivs):
            for i in range(1, len(ivs)):
                # check for adjacency/non-overlap
                # This constraint, combined with the start_i <= end_i constraint
                # of Intervals, guarantees that intervals are ordered
                assert ivs[i - 1].end_i == ivs[i].start_i, "intervals must be adjacent"

        self.interval_pairs = interval_pairs
        self.ref_intervals = [ip.ref for ip in self.interval_pairs]
        self.tgt_intervals = [ip.tgt for ip in self.interval_pairs]
        _validate_intervals(self.ref_intervals)
        _validate_intervals(self.tgt_intervals)
        self.ref_len = sum([iv.len for iv in self.ref_intervals])
        self.tgt_len = sum([iv.len for iv in self.tgt_intervals])

    @staticmethod
    def from_cigar(cigar):
        """
        :param cigar: a Compact Idiosyncratic Gapped Alignment Report string
        :type cigar: str.
        :returns: an IntervalMapper instance from the CIGAR string
        """
        return IntervalMapper(cigar_to_intervalpairs(cigar))

    def map_ref_to_tgt(self, start_i, end_i, max_extent=False):
        return self._map(self.ref_intervals, self.tgt_intervals, start_i, end_i, max_extent)

    def map_tgt_to_ref(self, start_i, end_i, max_extent=False):
        return self._map(self.tgt_intervals, self.ref_intervals, start_i, end_i, max_extent)

    @staticmethod
    def _map(from_ivs, to_ivs, from_start_i, from_end_i, max_extent):
        def iv_map(from_ivs, to_ivs, from_start_i, from_end_i, max_extent):
            """returns the <start,end> intervals indexes in which from_start_i and from_end_i occur"""
            # first look for 0-width interval that matches
            seil = [
                i for i, iv in enumerate(from_ivs)
                if iv.start_i == from_start_i and iv.end_i == from_end_i
            ]
            if len(seil) > 0:
                si = ei = seil[0]
            else:
                sil = [i for i, iv in enumerate(from_ivs) if iv.start_i <= from_start_i <= iv.end_i]
                eil = [i for i, iv in enumerate(from_ivs) if iv.start_i <= from_end_i <= iv.end_i]
                if len(sil) == 0 or len(eil) == 0:
                    raise HGVSInvalidIntervalError(
                        "start or end or both are beyond the bounds of transcript record")
                si, ei = (sil[0], eil[-1]) if max_extent else (sil[-1], eil[0])
            return si, ei

        def clip_to_iv(iv, pos):
            return max(iv.start_i, min(iv.end_i, pos))

        assert from_start_i <= from_end_i, "expected from_start_i <= from_end_i"
        try:
            si, ei = iv_map(from_ivs, to_ivs, from_start_i, from_end_i, max_extent)
        except ValueError:
            raise HGVSInvalidIntervalError("start_i,end_i interval out of bounds")
        to_start_i = clip_to_iv(to_ivs[si],
                                to_ivs[si].start_i + (from_start_i - from_ivs[si].start_i))
        to_end_i = clip_to_iv(to_ivs[ei], to_ivs[ei].end_i - (from_ivs[ei].end_i - from_end_i))
        return to_start_i, to_end_i


class CIGARElement(object):
    """represents elements of a CIGAR string and provides methods for
    determining the number of ref and tgt bases consumed by the
    operation"""

    __slots__ = ("len", "op")

    def __init__(self, len, op):
        self.len = len
        self.op = op

    @property
    def ref_len(self):
        """returns number of nt/aa consumed in reference sequence for this edit"""
        return self.len if self.op in "=INX" else 0

    @property
    def tgt_len(self):
        """returns number of nt/aa consumed in target sequence for this edit"""
        return self.len if self.op in "=DX" else 0


def cigar_to_intervalpairs(cigar):
    """For a given CIGAR string, return a list of (Interval,Interval)
    pairs.  The length of the returned list will be equal to the
    number of CIGAR operations
    """

    cigar_elem_re = re.compile(r"(?P<len>\d+)(?P<op>[=DIMNX])")
    ces = [
        CIGARElement(op=md["op"], len=int(md["len"]))
        for md in [m.groupdict() for m in cigar_elem_re.finditer(cigar)]
    ]
    ips = [None] * len(ces)
    ref_pos = tgt_pos = 0
    for i, ce in enumerate(ces):
        ips[i] = IntervalPair(
            Interval(ref_pos, ref_pos + ce.ref_len), Interval(tgt_pos, tgt_pos + ce.tgt_len))
        ref_pos += ce.ref_len
        tgt_pos += ce.tgt_len
    return ips


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
