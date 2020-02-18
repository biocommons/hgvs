# -*- coding: utf-8 -*-
"""Mapping positions between pairs of sequence alignments

The AlignmentMapper class is at the heart of mapping between aligned sequences.

"""

# Implementation note: HGVS doesn't have a 0. Counting is -3, -2, -1,
# 1, 2, 3 :-/ Coordinate calculations must take this discontinuity in
# c. positions into account.  The implementation of imaginary
# transcript positions creates a second discontinuity. (By analogy
# with c., n.0 is declared to not exist.)  The strategy used in this
# code is to use internal c0 and n0 coordinates, which include 0, for
# coordinate calculations and to translate these to c. and
# n. positions as needed.
#
#                imag.                                                 imag.
#              upstream     5' UTR           CDS          3' UTR      downstr
#                                     |>
#            - - - - - - ———————————— ||||||||||||||||| ——————————— - - - - - -
#                           a     b     C     D     E     f     g     h     i
#    c.        -4    -3    -2    -1  !  1     2     3  ! *1    *2    *3    *4
#    c0        -3    -2    -1     0     1     2     3     4     5     6     7
#    n0        -1     0     1     2     3     4     5     6     7     8     9
#    n.        -2    -1  !  1     2     3     4     5     6     7     8     9
#    g.   ... 123   124   125   126   127   128   129   130   131   132   133 ...
#



from __future__ import absolute_import, division, print_function, unicode_literals

from six.moves import range
from bioutils.coordinates import strand_int_to_pm

import hgvs.location

from hgvs import global_config
from hgvs.exceptions import HGVSError, HGVSUsageError, HGVSDataNotAvailableError, HGVSInvalidIntervalError
from hgvs.utils import build_tx_cigar
from hgvs.utils.cigarmapper import CIGARMapper
from hgvs.enums import Datum


class AlignmentMapper(object):
    """Maps hgvs location objects between genomic (g), non-coding (n) and
    cds (c) coordinates according to a CIGAR string.

    :param hdp: HGVS Data Provider Interface-compliant instance (see :class:`hgvs.dataproviders.interface.Interface`)
    :param str tx_ac: string representing transcript accession (e.g., NM_000551.2)
    :param str alt_ac: string representing the reference sequence accession (e.g., NC_000019.10)
    :param str alt_aln_method: string representing the alignment method; valid values depend on data source

    """
    __slots__ = ("tx_ac", "alt_ac", "alt_aln_method", "strand", "gc_offset", "cds_start_i",
                 "cds_end_i", "tgt_len", "cigarmapper", "ref_pos", "tgt_pos", "cigar_op")

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
                    "No transcript info".format(self=self))

            tx_exons = hdp.get_tx_exons(self.tx_ac, self.alt_ac, self.alt_aln_method)
            if tx_exons is None:
                raise HGVSDataNotAvailableError(
                    "AlignmentMapper(tx_ac={self.tx_ac}, "
                    "alt_ac={self.alt_ac}, alt_aln_method={self.alt_aln_method}): "
                    "No transcript exons".format(self=self))

            # hgvs-386: An assumption when building the cigar string
            # is that exons are adjacent. Assert that here.
            sorted_tx_exons = sorted(tx_exons, key=lambda e: e["ord"])
            for i in range(1, len(sorted_tx_exons)):
                if sorted_tx_exons[i - 1]["tx_end_i"] != sorted_tx_exons[i]["tx_start_i"]:
                    raise HGVSDataNotAvailableError(
                        "AlignmentMapper(tx_ac={self.tx_ac}, "
                        "alt_ac={self.alt_ac}, alt_aln_method={self.alt_aln_method}): "
                        "Exons {a} and {b} are not adjacent".format(self=self, a=i, b=i + 1))

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
                    "No transcript info".format(self=self))
            self.cds_start_i = tx_identity_info["cds_start_i"]
            self.cds_end_i = tx_identity_info["cds_end_i"]
            self.tgt_len = sum(tx_identity_info["lengths"])
            self.cigarmapper = None

        assert not (
            (self.cds_start_i is None) ^
            (self.cds_end_i is None)), "CDS start and end must both be defined or neither defined"


    def __str__(self):
        return "{self.__class__.__name__}: {self.tx_ac} ~ {self.alt_ac} ~ {self.alt_aln_method}; " \
               "{strand_pm} strand; offset={self.gc_offset}".format(
                    self=self, strand_pm=strand_int_to_pm(self.strand))

    def g_to_n(self, g_interval, strict_bounds=None):
        """convert a genomic (g.) interval to a transcript cDNA (n.) interval"""

        if strict_bounds is None:
            strict_bounds = global_config.mapping.strict_bounds

        grs, gre = g_interval.start.base - 1 - self.gc_offset, g_interval.end.base - 1 - self.gc_offset
        # frs, fre = (f)orward (r)na (s)tart & (e)nd; forward w.r.t. genome
        frs, frs_offset, frs_cigar = self.cigarmapper.map_ref_to_tgt(pos=grs, end="start", strict_bounds=strict_bounds)
        fre, fre_offset, fre_cigar = self.cigarmapper.map_ref_to_tgt(pos=gre, end="end", strict_bounds=strict_bounds)

        if self.strand == -1:
            frs, fre = self.tgt_len - fre - 1, self.tgt_len - frs - 1
            frs_offset, fre_offset = -fre_offset, -frs_offset
            
        # The returned interval would be uncertain when locating at alignment gaps
        return hgvs.location.BaseOffsetInterval(
            start=hgvs.location.BaseOffsetPosition(
                base=frs + 1, offset=frs_offset, datum=Datum.SEQ_START),
            end=hgvs.location.BaseOffsetPosition(
                base=fre + 1, offset=fre_offset, datum=Datum.SEQ_START),
            uncertain=frs_cigar in 'DI' or fre_cigar in 'DI')


    def n_to_g(self, n_interval, strict_bounds=None):
        """convert a transcript (n.) interval to a genomic (g.) interval"""

        if strict_bounds is None:
            strict_bounds = global_config.mapping.strict_bounds

        frs, fre = n_interval.start.base - 1, n_interval.end.base - 1
        start_offset, end_offset = n_interval.start.offset, n_interval.end.offset
        
        if self.strand == -1:
            fre, frs = self.tgt_len - frs - 1, self.tgt_len - fre - 1
            start_offset, end_offset = -end_offset, -start_offset

        # returns the genomic range start (grs) and end (gre)
        grs, _, grs_cigar = self.cigarmapper.map_tgt_to_ref(pos=frs, end="start", strict_bounds=strict_bounds)
        gre, _, gre_cigar = self.cigarmapper.map_tgt_to_ref(pos=fre, end="end", strict_bounds=strict_bounds)
        grs, gre = grs + self.gc_offset + 1, gre + self.gc_offset + 1
        gs, ge = grs + start_offset, gre + end_offset

        # The returned interval would be uncertain when locating at alignment gaps
        return hgvs.location.Interval(
            start=hgvs.location.SimplePosition(gs, uncertain=n_interval.start.uncertain),
            end=hgvs.location.SimplePosition(ge, uncertain=n_interval.end.uncertain),
            uncertain=grs_cigar in 'DI' or gre_cigar in 'DI')


    def n_to_c(self, n_interval, strict_bounds=None):
        """convert a transcript cDNA (n.) interval to a transcript CDS (c.) interval"""

        if strict_bounds is None:
            strict_bounds = global_config.mapping.strict_bounds

        if self.cds_start_i is None:    # cds_start_i defined iff cds_end_i defined; see assertion above
            raise HGVSUsageError(
                "CDS is undefined for {self.tx_ac}; cannot map to c. coordinate (non-coding transcript?)"
                .format(self=self))

        if strict_bounds and (n_interval.start.base <= 0 or n_interval.end.base > self.tgt_len):
            raise HGVSInvalidIntervalError(
                "The given coordinate is outside the bounds of the reference sequence.")


        # start
        if n_interval.start.base <= self.cds_start_i:
            cs = n_interval.start.base - (self.cds_start_i + 1)
            cs_datum = Datum.CDS_START
        elif n_interval.start.base > self.cds_start_i and n_interval.start.base <= self.cds_end_i:
            cs = n_interval.start.base - self.cds_start_i
            cs_datum = Datum.CDS_START
        else:
            cs = n_interval.start.base - self.cds_end_i
            cs_datum = Datum.CDS_END
        # end
        if n_interval.end.base <= self.cds_start_i:
            ce = n_interval.end.base - (self.cds_start_i + 1)
            ce_datum = Datum.CDS_START
        elif n_interval.end.base > self.cds_start_i and n_interval.end.base <= self.cds_end_i:
            ce = n_interval.end.base - self.cds_start_i
            ce_datum = Datum.CDS_START
        else:
            ce = n_interval.end.base - self.cds_end_i
            ce_datum = Datum.CDS_END

        c_interval = hgvs.location.BaseOffsetInterval(
            start=hgvs.location.BaseOffsetPosition(
                base=cs, offset=n_interval.start.offset, datum=cs_datum),
            end=hgvs.location.BaseOffsetPosition(
                base=ce, offset=n_interval.end.offset, datum=ce_datum),
            uncertain=n_interval.uncertain)
        return c_interval


    def c_to_n(self, c_interval, strict_bounds=None):
        """convert a transcript CDS (c.) interval to a transcript cDNA (n.) interval"""

        if strict_bounds is None:
            strict_bounds = global_config.mapping.strict_bounds

        if self.cds_start_i is None:    # cds_start_i defined iff cds_end_i defined; see assertion above
            raise HGVSUsageError(
                "CDS is undefined for {self.tx_ac}; cannot map from c. coordinate (non-coding transcript?)"
                .format(self=self))

        # start
        if c_interval.start.datum == Datum.CDS_START and c_interval.start.base < 0:
            n_start = c_interval.start.base + self.cds_start_i + 1
        elif c_interval.start.datum == Datum.CDS_START and c_interval.start.base > 0:
            n_start = c_interval.start.base + self.cds_start_i
        elif c_interval.start.datum == Datum.CDS_END:
            n_start = c_interval.start.base + self.cds_end_i
        # end
        if c_interval.end.datum == Datum.CDS_START and c_interval.end.base < 0:
            n_end = c_interval.end.base + self.cds_start_i + 1
        elif c_interval.end.datum == Datum.CDS_START and c_interval.end.base > 0:
            n_end = c_interval.end.base + self.cds_start_i
        elif c_interval.end.datum == Datum.CDS_END:
            n_end = c_interval.end.base + self.cds_end_i

        if strict_bounds and (n_start <= 0 or n_end > self.tgt_len):
            raise HGVSInvalidIntervalError(
                "The given coordinate is outside the bounds of the reference sequence.")
        
        n_interval = hgvs.location.BaseOffsetInterval(
            start=hgvs.location.BaseOffsetPosition(
                base=n_start, offset=c_interval.start.offset, datum=Datum.SEQ_START),
            end=hgvs.location.BaseOffsetPosition(
                base=n_end, offset=c_interval.end.offset, datum=Datum.SEQ_START),
            uncertain=c_interval.uncertain)
        return n_interval

    def g_to_c(self, g_interval, strict_bounds=None):
        """convert a genomic (g.) interval to a transcript CDS (c.) interval"""
        return self.n_to_c(self.g_to_n(g_interval), strict_bounds=strict_bounds)

    def c_to_g(self, c_interval, strict_bounds=None):
        """convert a transcript CDS (c.) interval to a genomic (g.) interval"""
        return self.n_to_g(self.c_to_n(c_interval), strict_bounds=strict_bounds)

    @property
    def is_coding_transcript(self):
        if ((self.cds_start_i is not None) ^ (self.cds_end_i is not None)):
            raise HGVSError("{self.tx_ac}: CDS start_i and end_i"
                            " must be both defined or both undefined".format(self=self))
        return self.cds_start_i is not None


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
