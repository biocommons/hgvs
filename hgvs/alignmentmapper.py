# -*- coding: utf-8 -*-
"""Mapping positions between pairs of sequence alignments

The AlignmentMapper class is at the heart of mapping between aligned sequences.

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import re
from six.moves import range
from bioutils.coordinates import strand_int_to_pm

import hgvs.location

from hgvs.exceptions import HGVSError, HGVSUsageError, HGVSDataNotAvailableError, HGVSInvalidIntervalError
from hgvs.utils import build_tx_cigar
from hgvs.enums import Datum

cigar_re = re.compile(r"(?P<len>\d+)(?P<op>[=DIMNX])")


class AlignmentMapper(object):
    """Provides coordinate (not variant) mapping operations between
    genomic (g), non-coding (n) and cds (c) coordinates according to a CIGAR.

    :param hdp: HGVS Data Provider Interface-compliant instance (see :class:`hgvs.dataproviders.interface.Interface`)
    :param str tx_ac: string representing transcript accession (e.g., NM_000551.2)
    :param str alt_ac: string representing the reference sequence accession (e.g., NC_000019.10)
    :param str alt_aln_method: string representing the alignment method; valid values depend on data source

    """
    __slots__ = ("tx_ac", "alt_ac", "alt_aln_method", "strand", "gc_offset", "cds_start_i",
                 "cds_end_i", "tgt_len", "cigar", "ref_pos", "tgt_pos", "cigar_op")

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
            self.cigar = build_tx_cigar(tx_exons, self.strand)
            self.ref_pos, self.tgt_pos, self.cigar_op = self._parse_cigar(self.cigar)
            self.tgt_len = self.tgt_pos[-1]
        else:
            # this covers the identity cases n <-> c
            tx_identity_info = hdp.get_tx_identity_info(self.tx_ac)
            if tx_identity_info is None:
                raise HGVSDataNotAvailableError(
                    "AlignmentMapper(tx_ac={self.tx_ac}, "
                    "alt_ac={self.alt_ac}, alt_aln_method={self.alt_aln_method}): "
                    "No transcript identity info".format(self=self))
            self.cds_start_i = tx_identity_info["cds_start_i"]
            self.cds_end_i = tx_identity_info["cds_end_i"]
            self.tgt_len = sum(tx_identity_info["lengths"])

        assert not (
            (self.cds_start_i is None) ^
            (self.cds_end_i is None)), "CDS start and end must both be defined or neither defined"

    def __str__(self):
        return "{self.__class__.__name__}: {self.tx_ac} ~ {self.alt_ac} ~ {self.alt_aln_method}; " \
               "{strand_pm} strand; offset={self.gc_offset}".format(
                    self=self, strand_pm=strand_int_to_pm(self.strand))

    def _parse_cigar(self, cigar):
        """For a given CIGAR string, return the start positions of
        each aligned segment in ref and tgt, and a list of CIGAR operators.
        """
        ces = [m.groupdict() for m in cigar_re.finditer(cigar)]
        ref_pos = [None] * len(ces)
        tgt_pos = [None] * len(ces)
        cigar_op = [None] * len(ces)
        ref_cur = tgt_cur = 0
        for i, ce in enumerate(ces):
            ref_pos[i] = ref_cur
            tgt_pos[i] = tgt_cur
            cigar_op[i] = ce["op"]
            step = int(ce["len"])
            if ce["op"] in "=MINX":
                ref_cur += step
            if ce["op"] in "=MDX":
                tgt_cur += step
        ref_pos.append(ref_cur)
        tgt_pos.append(tgt_cur)
        return ref_pos, tgt_pos, cigar_op

    def _map(self, from_pos, to_pos, pos, base):
        """Map position between aligned sequences

        Positions in this function are 0-based.
        """
        pos_i = -1
        while pos_i < len(self.cigar_op) and pos >= from_pos[pos_i + 1]:
            pos_i += 1

        if pos_i == -1 or pos_i == len(self.cigar_op):
            raise HGVSInvalidIntervalError("Position is beyond the bounds of transcript record")

        if self.cigar_op[pos_i] in "=MX":
            mapped_pos = to_pos[pos_i] + (pos - from_pos[pos_i])
            mapped_pos_offset = 0
        elif self.cigar_op[pos_i] in "DI":
            if base == "start":
                mapped_pos = to_pos[pos_i] - 1
            elif base == "end":
                mapped_pos = to_pos[pos_i]
            mapped_pos_offset = 0
        elif self.cigar_op[pos_i] == "N":
            if pos - from_pos[pos_i] + 1 <= from_pos[pos_i + 1] - pos:
                mapped_pos = to_pos[pos_i] - 1
                mapped_pos_offset = pos - from_pos[pos_i] + 1
            else:
                mapped_pos = to_pos[pos_i]
                mapped_pos_offset = -(from_pos[pos_i + 1] - pos)

        return mapped_pos, mapped_pos_offset, self.cigar_op[pos_i]

    def g_to_n(self, g_interval):
        """convert a genomic (g.) interval to a transcript cDNA (n.) interval"""

        grs, gre = g_interval.start.base - 1 - self.gc_offset, g_interval.end.base - 1 - self.gc_offset
        # frs, fre = (f)orward (r)na (s)tart & (e)nd; forward w.r.t. genome
        frs, frs_offset, frs_cigar = self._map(
            from_pos=self.ref_pos, to_pos=self.tgt_pos, pos=grs, base="start")
        fre, fre_offset, fre_cigar = self._map(
            from_pos=self.ref_pos, to_pos=self.tgt_pos, pos=gre, base="end")

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

    def n_to_g(self, n_interval):
        """convert a transcript (n.) interval to a genomic (g.) interval"""

        frs, fre = n_interval.start.base - 1, n_interval.end.base - 1
        start_offset, end_offset = n_interval.start.offset, n_interval.end.offset

        if self.strand == -1:
            fre, frs = self.tgt_len - frs - 1, self.tgt_len - fre - 1
            start_offset, end_offset = -end_offset, -start_offset

        # returns the genomic range start (grs) and end (gre)
        grs, _, grs_cigar = self._map(
            from_pos=self.tgt_pos, to_pos=self.ref_pos, pos=frs, base="start")
        gre, _, gre_cigar = self._map(
            from_pos=self.tgt_pos, to_pos=self.ref_pos, pos=fre, base="end")
        grs, gre = grs + self.gc_offset + 1, gre + self.gc_offset + 1
        gs, ge = grs + start_offset, gre + end_offset

        # The returned interval would be uncertain when locating at alignment gaps
        return hgvs.location.Interval(
            start=hgvs.location.SimplePosition(gs, uncertain=n_interval.start.uncertain),
            end=hgvs.location.SimplePosition(ge, uncertain=n_interval.end.uncertain),
            uncertain=grs_cigar in 'DI' or gre_cigar in 'DI')

    def n_to_c(self, n_interval):
        """convert a transcript cDNA (n.) interval to a transcript CDS (c.) interval"""

        if self.cds_start_i is None:    # cds_start_i defined iff cds_end_i defined; see assertion above
            raise HGVSUsageError(
                "CDS is undefined for {self.tx_ac}; cannot map to c. coordinate (non-coding transcript?)"
                .format(self=self))
        if n_interval.start.base <= 0 or n_interval.end.base > self.tgt_len:
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

    def c_to_n(self, c_interval):
        """convert a transcript CDS (c.) interval to a transcript cDNA (n.) interval"""

        if self.cds_start_i is None:    # cds_start_i defined iff cds_end_i defined; see assertion above
            raise HGVSUsageError(
                "CDS is undefined for {self.tx_ac}; cannot map from c. coordinate (non-coding transcript?)"
                .format(self=self))

        # start
        if c_interval.start.datum == Datum.CDS_START and c_interval.start.base < 0:
            r_start = c_interval.start.base + self.cds_start_i + 1
        elif c_interval.start.datum == Datum.CDS_START and c_interval.start.base > 0:
            r_start = c_interval.start.base + self.cds_start_i
        elif c_interval.start.datum == Datum.CDS_END:
            r_start = c_interval.start.base + self.cds_end_i
        # end
        if c_interval.end.datum == Datum.CDS_START and c_interval.end.base < 0:
            r_end = c_interval.end.base + self.cds_start_i + 1
        elif c_interval.end.datum == Datum.CDS_START and c_interval.end.base > 0:
            r_end = c_interval.end.base + self.cds_start_i
        elif c_interval.end.datum == Datum.CDS_END:
            r_end = c_interval.end.base + self.cds_end_i

        if r_start <= 0 or r_end > self.tgt_len:
            raise HGVSInvalidIntervalError(
                "The given coordinate is outside the bounds of the reference sequence.")

        n_interval = hgvs.location.BaseOffsetInterval(
            start=hgvs.location.BaseOffsetPosition(
                base=r_start, offset=c_interval.start.offset, datum=Datum.SEQ_START),
            end=hgvs.location.BaseOffsetPosition(
                base=r_end, offset=c_interval.end.offset, datum=Datum.SEQ_START),
            uncertain=c_interval.uncertain)
        return n_interval

    def g_to_c(self, g_interval):
        """convert a genomic (g.) interval to a transcript CDS (c.) interval"""
        return self.n_to_c(self.g_to_n(g_interval))

    def c_to_g(self, c_interval):
        """convert a transcript CDS (c.) interval to a genomic (g.) interval"""
        return self.n_to_g(self.c_to_n(c_interval))

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
