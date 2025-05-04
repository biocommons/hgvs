# -*- coding: utf-8 -*-
"""Mapping positions between pairs of sequence alignments

The AlignmentMapper class is at the heart of mapping between aligned sequences.

"""

# Implementation note re: "no-zero correction": HGVS doesn't have a
# 0. Counting is -3, -2, -1, 1, 2, 3 :-/ Coordinate calculations must
# take this discontinuity in c. positions into account.  The
# implementation of imaginary transcript positions creates a second
# discontinuity. (By analogy with c., n.0 is declared to not exist.)
# The strategy used in this code is to use internal c0 and n0
# coordinates, which include 0, for coordinate calculations and to
# translate these to c. and n. positions as needed.
#
#                imag.                                                 imag.
#              upstream     5' UTR           CDS          3' UTR      downstr
#                                     |>
#            - - - - - - ———————————— ||||||||||||||||| ——————————— - - - - - -
#                           a     b     C     D     E     f     g     h     i
#    c.        -4    -3    -2    -1  !  1     2     3  ! *1    *2    *3    *4
#    c0        -4    -3    -2    -1     0     1     2     3     4     5     6
#    n0        -2    -1     0     1     2     3     4     5     6     7     8
#    n.        -2    -1  !  1     2     3     4     5     6     7     8     9
#    g.   ... 123   124   125   126   127   128   129   130   131   132   133 ...
#

from __future__ import absolute_import, division, print_function, unicode_literals

from typing import Optional

from bioutils.coordinates import strand_int_to_pm

import hgvs.location
from hgvs import global_config
from hgvs.enums import Datum
from hgvs.exceptions import (
    HGVSDataNotAvailableError,
    HGVSError,
    HGVSInvalidIntervalError,
    HGVSUsageError,
)
from hgvs.location import Interval, BaseOffsetInterval
from hgvs.utils import build_tx_cigar
from hgvs.utils.cigarmapper import CIGARMapper


def _zbc_to_hgvs(i: int):
    """Convert zero-based coordinate to hgvs (1 based, missing zero)"""
    if i >= 0:
        i += 1
    return i


def _hgvs_to_zbc(i: int):
    """Convert hgvs (1 based, missing zero)"""
    if i >= 1:
        i -= 1
    return i


class AlignmentMapper:
    """Maps hgvs location objects between genomic (g), non-coding (n) and
    cds (c) coordinates according to a CIGAR string.

    :param hdp: HGVS Data Provider Interface-compliant instance (see :class:`hgvs.dataproviders.interface.Interface`)
    :param str tx_ac: string representing transcript accession (e.g., NM_000551.2)
    :param str alt_ac: string representing the reference sequence accession (e.g., NC_000019.10)
    :param str alt_aln_method: string representing the alignment method; valid values depend on data source

    """

    __slots__ = (
        "tx_ac",
        "alt_ac",
        "alt_aln_method",
        "strand",
        "gc_offset",
        "cds_start_i",
        "cds_end_i",
        "tgt_len",
        "cigarmapper",
        "ref_pos",
        "tgt_pos",
        "cigar_op",
    )

    def __init__(self, hdp, tx_ac, alt_ac, alt_aln_method):
        self.tx_ac = tx_ac
        self.alt_ac = alt_ac
        self.alt_aln_method = alt_aln_method

        if self.alt_aln_method != "transcript":
            tx_info = hdp.get_tx_info(self.tx_ac, self.alt_ac, self.alt_aln_method)
            if tx_info is None:
                raise HGVSDataNotAvailableError(
                    "AlignmentMapper(tx_ac={self.tx_ac}, "
                    "alt_ac={self.alt_ac}, alt_aln_method={self.alt_aln_method}): "
                    "No transcript info".format(self=self)
                )

            tx_exons = hdp.get_tx_exons(self.tx_ac, self.alt_ac, self.alt_aln_method)
            if tx_exons is None:
                raise HGVSDataNotAvailableError(
                    "AlignmentMapper(tx_ac={self.tx_ac}, "
                    "alt_ac={self.alt_ac}, alt_aln_method={self.alt_aln_method}): "
                    "No transcript exons".format(self=self)
                )

            # hgvs-386: An assumption when building the cigar string
            # is that exons are adjacent. Assert that here.
            sorted_tx_exons = sorted(tx_exons, key=lambda e: e["ord"])
            for i in range(1, len(sorted_tx_exons)):
                if (
                    sorted_tx_exons[i - 1]["tx_end_i"]
                    != sorted_tx_exons[i]["tx_start_i"]
                ):
                    raise HGVSDataNotAvailableError(
                        "AlignmentMapper(tx_ac={self.tx_ac}, "
                        "alt_ac={self.alt_ac}, alt_aln_method={self.alt_aln_method}): "
                        "Exons {a} and {b} are not adjacent".format(
                            self=self, a=i, b=i + 1
                        )
                    )

            self.strand = tx_exons[0]["alt_strand"]
            self.gc_offset = tx_exons[0]["alt_start_i"]
            self.cds_start_i = tx_info["cds_start_i"]
            self.cds_end_i = tx_info["cds_end_i"]

            cigar = build_tx_cigar(tx_exons, self.strand)
            self.cigarmapper = CIGARMapper(cigar)
            self.tgt_len = self.cigarmapper.tgt_len

        else:
            # this covers the identity cases n <-> c
            tx_identity_info = hdp.get_tx_identity_info(self.tx_ac)
            if tx_identity_info is None:
                raise HGVSDataNotAvailableError(
                    "AlignmentMapper(tx_ac={self.tx_ac}, "
                    "alt_ac={self.alt_ac}, alt_aln_method={self.alt_aln_method}): "
                    "No transcript info".format(self=self)
                )
            self.cds_start_i = tx_identity_info["cds_start_i"]
            self.cds_end_i = tx_identity_info["cds_end_i"]
            self.tgt_len = sum(tx_identity_info["lengths"])
            self.cigarmapper = None

        assert not ((self.cds_start_i is None) ^ (self.cds_end_i is None)), (
            "CDS start and end must both be defined or neither defined"
        )

    def __str__(self):
        return (
            "{self.__class__.__name__}: {self.tx_ac} ~ {self.alt_ac} ~ {self.alt_aln_method}; "
            "{strand_pm} strand; offset={self.gc_offset}".format(
                self=self, strand_pm=strand_int_to_pm(self.strand)
            )
        )

    def _g_to_n_interval(
        self, g_interval: Interval, strict_bounds: Optional[bool] = None
    ) -> BaseOffsetInterval:
        """convert a genomic (g.) interval to a transcript cDNA (n.) interval"""

        ss, se = self._get_start_end(g_interval.start)
        es, ee = self._get_start_end(g_interval.end)

        print(f"g_to_n: g_interval {g_interval} ss {ss} se {se} es {es} ee {ee}")

        grs_left = ss.base - 1 - self.gc_offset
        grs_right = se.base - 1 - self.gc_offset
        gre_left = es.base - 1 - self.gc_offset
        gre_right = ee.base - 1 - self.gc_offset

        # frs, fre = (f)orward (r)na (s)tart & (e)nd; forward w.r.t. genome
        (
            forward_rna_start_left,
            forward_rna_start_offset_left,
            forward_rna_start_cigar_left,
        ) = self.cigarmapper.map_ref_to_tgt(
            pos=grs_left, end="start", strict_bounds=strict_bounds
        )
        (
            forward_rna_start_right,
            forward_rna_start_offset_right,
            forward_rna_start_cigar_right,
        ) = self.cigarmapper.map_ref_to_tgt(
            pos=grs_right, end="start", strict_bounds=strict_bounds
        )
        (
            forward_rna_end_left,
            forward_rna_end_offset_left,
            forward_rna_end_cigar_left,
        ) = self.cigarmapper.map_ref_to_tgt(
            pos=gre_left, end="end", strict_bounds=strict_bounds
        )
        (
            forward_rna_end_right,
            forward_rna_end_offset_right,
            forward_rna_end_cigar_right,
        ) = self.cigarmapper.map_ref_to_tgt(
            pos=gre_right, end="end", strict_bounds=strict_bounds
        )

        if self.strand == -1:
            forward_rna_start_left = self.tgt_len - 1 - forward_rna_end_right
            forward_rna_start_right = self.tgt_len - 1 - forward_rna_end_left
            forward_rna_end_left = self.tgt_len - 1 - forward_rna_start_right
            forward_rna_end_right = self.tgt_len - 1 - forward_rna_start_left

            forward_rna_start_offset_left = -forward_rna_end_offset_right
            forward_rna_start_offset_right = -forward_rna_end_offset_left
            forward_rna_end_offset_right = -forward_rna_start_offset_left
            forward_rna_end_offset_left = -forward_rna_start_offset_right

        # make this an interval, not a point
        start_left = hgvs.location.BaseOffsetPosition(
            base=_zbc_to_hgvs(forward_rna_start_left),
            offset=forward_rna_start_offset_left,
            datum=Datum.SEQ_START,
            uncertain=forward_rna_start_cigar_left in "DI",
        )
        start_right = hgvs.location.BaseOffsetPosition(
            base=_zbc_to_hgvs(forward_rna_start_right),
            offset=forward_rna_start_offset_right,
            datum=Datum.SEQ_START,
            uncertain=forward_rna_start_cigar_right in "DI",
        )

        start = hgvs.location.BaseOffsetInterval(
            start=start_left,
            end=start_right,
            uncertain=g_interval.start.uncertain,
        )
        end_left = hgvs.location.BaseOffsetPosition(
            base=_zbc_to_hgvs(forward_rna_end_left),
            offset=forward_rna_end_offset_left,
            datum=Datum.SEQ_START,
            uncertain=forward_rna_end_cigar_left in "DI",
        )
        end_right = hgvs.location.BaseOffsetPosition(
            base=_zbc_to_hgvs(forward_rna_end_right),
            offset=forward_rna_end_offset_right,
            datum=Datum.SEQ_START,
            uncertain=forward_rna_end_cigar_right in "DI",
        )
        end = hgvs.location.BaseOffsetInterval(
            start=end_left,
            end=end_right,
            uncertain=g_interval.end.uncertain,
        )

        # The returned interval would be uncertain when locating at alignment gaps
        # of if the initial interval was uncertain
        return hgvs.location.Interval(
            start=start,
            end=end,
            uncertain=g_interval.uncertain,
        )

    def g_to_n(
        self,
        g_interval: Interval,
        strict_bounds: Optional[bool] = None,
        imprecise_inner_interval_only: bool | None = None,
    ) -> BaseOffsetInterval:
        """convert a genomic (g.) interval to a transcript cDNA (n.) interval"""

        if strict_bounds is None:
            strict_bounds = global_config.mapping.strict_bounds

        if imprecise_inner_interval_only is None:
            imprecise_inner_interval_only = (
                global_config.g_to_c.imprecise_inner_interval_only
            )

        # add config parameter:
        print(
            f"g_to_n: g_interval {g_interval} imprecise_inner_interval_only {imprecise_inner_interval_only}"
        )
        if (
            not imprecise_inner_interval_only
            and isinstance(g_interval.start, Interval)
            and isinstance(g_interval.end, Interval)
        ):
            return self._g_to_n_interval(g_interval, strict_bounds)

        print(
            f"check ranges {self.strand} {g_interval.start} {g_interval.start.uncertain} {imprecise_inner_interval_only}"
        )
        # in case of uncertain ranges, we fall back to the inner (more confident) interval
        if imprecise_inner_interval_only and g_interval.start.uncertain:
            grs = g_interval.start.end.base - 1 - self.gc_offset
        elif isinstance(g_interval.start, Interval):
            grs = g_interval.start.end.base - 1 - self.gc_offset
        else:
            grs = g_interval.start.base - 1 - self.gc_offset

        print(f"g_to_n: grs {grs} {imprecise_inner_interval_only}")

        if imprecise_inner_interval_only and g_interval.end.uncertain:
            gre = g_interval.end.start.base - 1 - self.gc_offset
        elif isinstance(g_interval.end, Interval):
            gre = g_interval.end.start.base - 1 - self.gc_offset
        else:
            gre = g_interval.end.base - 1 - self.gc_offset

        print(f"g_to_n: gre {gre} {imprecise_inner_interval_only}")

        forward_rna_start, forward_rna_start_offset, forward_rna_start_cigar = (
            self.cigarmapper.map_ref_to_tgt(
                pos=grs, end="start", strict_bounds=strict_bounds
            )
        )
        forward_rna_end, forward_rna_end_offset, forward_rna_end_cigar = (
            self.cigarmapper.map_ref_to_tgt(
                pos=gre, end="end", strict_bounds=strict_bounds
            )
        )
        if self.strand == -1:
            forward_rna_start = self.tgt_len - 1 - forward_rna_end
            forward_rna_end = self.tgt_len - 1 - forward_rna_start

            forward_rna_start_offset = -forward_rna_end_offset
            forward_rna_end_offset = -forward_rna_start_offset
            print(
                f"g_to_n: strand -1 {forward_rna_start} {forward_rna_end} {forward_rna_start_offset} {forward_rna_end_offset}"
            )

        # The returned interval would be uncertain when locating at alignment gaps
        # of if the initial interval was uncertain
        return hgvs.location.BaseOffsetInterval(
            start=hgvs.location.BaseOffsetPosition(
                base=_zbc_to_hgvs(forward_rna_start),
                offset=forward_rna_start_offset,
                datum=Datum.SEQ_START,
                uncertain=g_interval.start.uncertain,
            ),
            end=hgvs.location.BaseOffsetPosition(
                base=_zbc_to_hgvs(forward_rna_end),
                offset=forward_rna_end_offset,
                datum=Datum.SEQ_START,
                uncertain=g_interval.end.uncertain,
            ),
            uncertain=forward_rna_start_cigar in "DI" or forward_rna_end_cigar in "DI",
        )

    def _get_start_end(self, var):
        print(f"get_start_end {var} {type(var)}")

        if isinstance(var, hgvs.location.Interval):
            print(
                f"get_start_end Interval {var.start} {var.end} {type(var.start)} {type(var.end)}"
            )
            _, s = self._get_start_end(var.start)
            e, _ = self._get_start_end(var.end)

            if not s.base:
                s = e
            if not e.base:
                e = s

            return s, e

        if isinstance(var, hgvs.location.BaseOffsetPosition) or isinstance(
            var, hgvs.location.SimplePosition
        ):
            s = var
            e = var

            return s, e

        if isinstance(var.posedit, hgvs.location.BaseOffsetInterval):
            s = var.posedit.start
            e = var.posedit.end

        if isinstance(var.posedit.pos, hgvs.location.Interval):
            s = var.posedit.pos.start.start
        else:
            s = var.posedit.pos.start
        if isinstance(var.posedit.pos.end, hgvs.location.Interval):
            e = var.posedit.pos.end.end
        else:
            e = var.posedit.pos.end

        return s, e

    def n_to_g(self, n_interval, strict_bounds=None) -> Interval:
        """convert a transcript (n.) interval to a genomic (g.) interval"""

        if strict_bounds is None:
            strict_bounds = global_config.mapping.strict_bounds

        s, e = self._get_start_end(n_interval)

        print(
            f"n_to_g: {n_interval} {s} {e} {type(n_interval)} {s.datum} {s.uncertain}"
        )
        frs = _hgvs_to_zbc(s.base)
        start_offset = s.offset
        fre = _hgvs_to_zbc(e.base)
        end_offset = e.offset

        if self.strand == -1:
            fre, frs = self.tgt_len - 1 - frs, self.tgt_len - 1 - fre
            start_offset, end_offset = -end_offset, -start_offset

        # returns the genomic range start (grs) and end (gre)
        grs, _, grs_cigar = self.cigarmapper.map_tgt_to_ref(
            pos=frs, end="start", strict_bounds=strict_bounds
        )
        gre, _, gre_cigar = self.cigarmapper.map_tgt_to_ref(
            pos=fre, end="end", strict_bounds=strict_bounds
        )
        grs, gre = grs + self.gc_offset + 1, gre + self.gc_offset + 1
        gs, ge = grs + start_offset, gre + end_offset

        if n_interval.start.uncertain:
            start = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(uncertain=False),
                end=hgvs.location.SimplePosition(gs, uncertain=False),
                uncertain=True,
            )
        else:
            start = hgvs.location.SimplePosition(gs, uncertain=False)

        if n_interval.end.uncertain:
            end = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(ge, uncertain=False),
                end=hgvs.location.SimplePosition(uncertain=False),
                uncertain=True,
            )
        else:
            end = hgvs.location.SimplePosition(ge, uncertain=False)

        # The returned interval would be uncertain when locating at alignment gaps
        return hgvs.location.Interval(
            start=start,
            end=end,
            uncertain=grs_cigar in "DI" or gre_cigar in "DI",
        )

    def n_to_c(
        self,
        n_interval: Interval,
        strict_bounds: bool | None = None,
        imprecise_inner_interval_only: bool | None = None,
    ):
        """convert a transcript cDNA (n.) interval to a transcript CDS (c.) interval"""

        if strict_bounds is None:
            strict_bounds = global_config.mapping.strict_bounds
        if imprecise_inner_interval_only is None:
            imprecise_inner_interval_only = (
                global_config.g_to_c.imprecise_inner_interval_only
            )

        if (
            self.cds_start_i is None
        ):  # cds_start_i defined iff cds_end_i defined; see assertion above
            raise HGVSUsageError(
                "CDS is undefined for {self.tx_ac}; cannot map to c. coordinate (non-coding transcript?)".format(
                    self=self
                )
            )

        if isinstance(n_interval, BaseOffsetInterval):
            start = n_interval.start.base
            end = n_interval.end.base
        elif isinstance(n_interval, Interval):
            start = n_interval.start.start.base
            end = n_interval.end.end.base

        if strict_bounds and (start <= 0 or end > self.tgt_len):
            raise HGVSInvalidIntervalError(
                "The given coordinate is outside the bounds of the reference sequence."
            )

        def pos_n_to_c(pos) -> hgvs.location.BaseOffsetPosition:
            if pos.base <= self.cds_start_i:
                c = pos.base - self.cds_start_i - (1 if pos.base > 0 else 0)
                c_datum = Datum.CDS_START
            elif pos.base > self.cds_start_i and pos.base <= self.cds_end_i:
                c = pos.base - self.cds_start_i
                c_datum = Datum.CDS_START
            else:
                c = pos.base - self.cds_end_i
                c_datum = Datum.CDS_END
            return hgvs.location.BaseOffsetPosition(
                base=c, offset=pos.offset, datum=c_datum, uncertain=pos.uncertain
            )

        def interval_n_to_c(interval: Interval) -> hgvs.location.BaseOffsetInterval:
            istart = interval.start
            c_start = pos_n_to_c(istart)

            iend = interval.end
            c_end = pos_n_to_c(iend)
            return hgvs.location.BaseOffsetInterval(
                start=c_start,
                end=c_end,
                uncertain=interval.uncertain,
            )

        baseoffset = True
        if isinstance(n_interval.start, Interval):
            c_start = interval_n_to_c(n_interval.start)
            baseoffset = False
        else:
            c_start = pos_n_to_c(n_interval.start)

        if isinstance(n_interval.end, Interval):
            c_end = interval_n_to_c(n_interval.end)
            baseoffset = False
        else:
            # pos
            c_end = pos_n_to_c(n_interval.end)

        if baseoffset:
            c_interval = hgvs.location.BaseOffsetInterval(
                start=c_start,
                end=c_end,
                uncertain=n_interval.uncertain,
            )
        else:
            c_interval = hgvs.location.Interval(
                start=c_start,
                end=c_end,
                uncertain=n_interval.uncertain,
            )
        return c_interval

    def c_to_n(self, c_interval, strict_bounds=None):
        """convert a transcript CDS (c.) interval to a transcript cDNA (n.) interval"""

        if strict_bounds is None:
            strict_bounds = global_config.mapping.strict_bounds

        if self.cds_start_i is None:
            raise HGVSUsageError(
                "CDS is undefined for {self.tx_ac}; this accession appears to be for a non-coding transcript".format(
                    self=self
                )
            )

        def pos_c_to_n(pos):
            if pos.datum == Datum.CDS_START:
                n = pos.base + self.cds_start_i
                if pos.base < 0:  # correct for lack of c.0 coordinate
                    n += 1
            elif pos.datum == Datum.CDS_END:
                n = pos.base + self.cds_end_i
            if n <= 0:  # correct for lack of n.0 coordinate
                n -= 1
            if n <= 0 or n > self.tgt_len:
                if strict_bounds:
                    raise HGVSInvalidIntervalError(
                        f"c.{pos} coordinate is out of bounds"
                    )

            return hgvs.location.BaseOffsetPosition(
                base=n,
                offset=pos.offset,
                datum=Datum.SEQ_START,
                uncertain=pos.uncertain,
            )

        def interval_c_to_n(interval: Interval) -> hgvs.location.BaseOffsetInterval:
            istart = interval.start
            n_start = pos_c_to_n(istart)

            iend = interval.end
            n_end = pos_c_to_n(iend)

            return hgvs.location.BaseOffsetInterval(
                start=n_start,
                end=n_end,
                uncertain=interval.uncertain,
            )

        print(f"c_to_n {type(c_interval)}")
        print(f"c_to_n {c_interval}")

        baseoffset = True
        if isinstance(c_interval.start, Interval):
            n_start = interval_c_to_n(c_interval.start)
            baseoffset = False
            print(
                f"c to n: I c_interval start {c_interval.start} c_interval end {c_interval.end}"
            )
            print(f"c to n: I start {n_start}")
        else:
            print(
                f"c to n: c_interval start {c_interval.start} c_interval end {c_interval.end}"
            )
            print(
                f"c to n: start {pos_c_to_n(c_interval.start)} end {pos_c_to_n(c_interval.end)}"
            )
            n_start = pos_c_to_n(c_interval.start)

        if isinstance(c_interval.end, Interval):
            n_end = interval_c_to_n(c_interval.end)
            baseoffset = False
        else:
            n_end = pos_c_to_n(c_interval.end)

        if baseoffset:
            n_interval = hgvs.location.BaseOffsetInterval(
                start=n_start,
                end=n_end,
                uncertain=c_interval.uncertain,
            )
        else:
            n_interval = hgvs.location.Interval(
                start=n_start,
                end=n_end,
                uncertain=c_interval.uncertain,
            )

        return n_interval

    def g_to_c(
        self,
        g_interval,
        strict_bounds: bool | None = None,
        imprecise_inner_interval_only: bool | None = None,
    ):
        """convert a genomic (g.) interval to a transcript CDS (c.) interval"""
        return self.n_to_c(
            self.g_to_n(
                g_interval, imprecise_inner_interval_only=imprecise_inner_interval_only
            ),
            strict_bounds=strict_bounds,
            imprecise_inner_interval_only=imprecise_inner_interval_only,
        )

    def c_to_g(self, c_interval, strict_bounds=None):
        """convert a transcript CDS (c.) interval to a genomic (g.) interval"""
        return self.n_to_g(self.c_to_n(c_interval), strict_bounds=strict_bounds)

    @property
    def is_coding_transcript(self):
        if (self.cds_start_i is not None) ^ (self.cds_end_i is not None):
            raise HGVSError(
                "{self.tx_ac}: CDS start_i and end_i"
                " must be both defined or both undefined".format(self=self)
            )
        return self.cds_start_i is not None

    def g_interval_is_inbounds(self, ival):
        grs = ival.start.base - 1 - self.gc_offset
        gre = ival.end.base - 1 - self.gc_offset
        return grs >= 0 and gre <= self.cigarmapper.ref_len


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
