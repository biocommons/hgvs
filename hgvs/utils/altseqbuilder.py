# -*- coding: utf-8 -*-
"""Utility to insert an hgvs variant into a transcript sequence.
Generates a record corresponding to the modified transcript sequence,
along with annotations for use in conversion to an hgvsp tag.
Used in hgvsc to hgvsp conversion.

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import math

from Bio.Seq import Seq
from bioutils.sequences import reverse_complement

from ..edit import (NARefAlt, Dup, Inv, Repeat)
from ..enums import Datum
import six

DBG = False

_logger = logging.getLogger(__name__)


class AltTranscriptData(object):
    def __init__(self,
                 seq,
                 cds_start,
                 cds_stop,
                 is_frameshift,
                 variant_start_aa,
                 accession,
                 is_substitution=False,
                 is_ambiguous=False):
        """Create a variant sequence using inputs from VariantInserter
        :param seq: DNA sequence wiith variant incorporated
        :type seq: str or list
        :param cds_start: coding sequence start (1-based)
        :type cds_start: int
        :param cds_stop: coding sequence stop (1-based)
        :type cds_stop: int
        :param protein_accession: protein accession, e.g. NP_999999.2
        :type protein_accession: str
        :param is_frameshift: is this variant a frameshift
        :type is_frameshift: bool
        :param variant_start_aa: AA start index (1-based) for this variant
        :type variant_start_aa: int
        :param is_substitution: flag if this is a substitution AA variant
        :type is_substitution: bool
        :param is_ambiguous: flag if variant is "?"
        :type is_ambiguous: bool
        :return variant sequence data
        :rtype attrs
        """

        if len(seq) > 0:
            if isinstance(seq, six.string_types):
                seq = list(seq)
            seq_cds = seq[cds_start - 1:]
            if len(seq_cds) % 3 != 0:    # padding so biopython won't complain during the conversion
                seq_cds.extend(['N'] * ((3 - len(seq_cds) % 3) % 3))
            seq_cds = ''.join(seq_cds)
            seq_aa = str(Seq(seq_cds).translate())
            stop_pos = seq_aa[:(cds_stop - cds_start + 1) // 3].rfind("*")
            if stop_pos == -1:
                stop_pos = seq_aa.find("*")
            if stop_pos != -1:
                seq_aa = seq_aa[:stop_pos + 1]
        else:
            seq_aa = []

        self.transcript_sequence = ''.join(seq)
        self.aa_sequence = seq_aa
        self.cds_start = cds_start
        self.cds_stop = cds_stop
        self.protein_accession = accession
        self.is_frameshift = is_frameshift
        self.variant_start_aa = variant_start_aa
        self.frameshift_start = None
        self.is_substitution = is_substitution
        self.is_ambiguous = is_ambiguous


class AltSeqBuilder(object):

    EXON = "exon"
    INTRON = "intron"
    F_UTR = "five utr"
    T_UTR = "three utr"
    WHOLE_GENE = "whole gene"

    def __init__(self, var_c, transcript_data):
        """Constructor

        :param var_c: representation of hgvs variant
        :type var_c: SequenceVariant
        :param transcript_data: representation of transcript
        :type transcript_data: attrs

        """
        self._var_c = var_c
        self._transcript_data = transcript_data
        if DBG:
            print(transcript_data.transcript_sequence)

        # check reference for special characteristics
        self._ref_has_multiple_stops = self._transcript_data.aa_sequence.count("*") > 1

    def build_altseq(self):
        """given a variant and a sequence, incorporate the variant and return the new sequence

        Data structure returned is analogous to the data structure used to return the variant sequence,
        but with an additional parameter denoting the start of a frameshift that should affect all bases
        downstream.

        :returns variant sequence data
        :rtype list of dictionaries
        """
        NOT_CDS = "not_cds_variant"
        WHOLE_GENE_DELETED = "whole_gene_deleted"

        type_map = {
            NARefAlt: self._incorporate_delins,
            Dup: self._incorporate_dup,
            Inv: self._incorporate_inv,
            Repeat: self._incorporate_repeat,
            NOT_CDS: self._create_alt_equals_ref_noncds,
            WHOLE_GENE_DELETED: self._create_no_protein
        }

        # should loop over each allele rather than assume only 1 variant; return a list for now
        alt_data = []

        variant_location = self._get_variant_region()

        if variant_location == self.EXON:
            edit_type = type(self._var_c.posedit.edit)
        elif variant_location == self.INTRON:
            edit_type = NOT_CDS
        elif variant_location == self.T_UTR:
            edit_type = NOT_CDS
        elif variant_location == self.F_UTR:
            # TODO: handle case where variant introduces a Met (new start)
            edit_type = NOT_CDS
        elif variant_location == self.WHOLE_GENE:
            if self._var_c.posedit.edit.type == "del":
                edit_type = WHOLE_GENE_DELETED
            elif self._var_c.posedit.edit.type == "dup":
                _logger.warning(
                    "Whole-gene duplication; consequence assumed to not affect protein product")
                edit_type = NOT_CDS
            elif self._var_c.posedit.edit.type == "inv":
                _logger.warning(
                    "Whole-gene inversion; consequence assumed to not affect protein product")
                edit_type = NOT_CDS
            else:
                edit_type = NOT_CDS
        else:    # should never get here
            raise ValueError("value_location = {}".format(variant_location))

        try:
            this_alt_data = type_map[edit_type]()
        except KeyError:
            raise NotImplementedError("c to p translation unsupported for {} type {}".format(
                self._var_c, edit_type))

        # get the start of the "terminal" frameshift (i.e. one never "cancelled out")
        this_alt_data = self._get_frameshift_start(this_alt_data)
        alt_data.append(this_alt_data)
        if DBG:
            print(this_alt_data.transcript_sequence)

        return alt_data

    def _get_variant_region(self):
        """Categorize variant by location in transcript (5'utr, exon, intron, 3'utr)

        :return "exon", "intron", "five_utr", "three_utr", "whole_gene"
        :rtype str
        """
        if self._var_c.posedit.pos.start.datum == Datum.CDS_END and self._var_c.posedit.pos.end.datum == Datum.CDS_END:
            result = self.T_UTR
        elif self._var_c.posedit.pos.start.base < 0 and self._var_c.posedit.pos.end.base < 0:
            result = self.F_UTR
        elif self._var_c.posedit.pos.start.base < 0 and self._var_c.posedit.pos.end.datum == Datum.CDS_END:
            result = self.WHOLE_GENE
        elif self._var_c.posedit.pos.start.offset != 0 or self._var_c.posedit.pos.end.offset != 0:
            # leave out anything intronic for now
            result = self.INTRON
        else:    # anything else that contains an exon
            result = self.EXON
        return result

    # def _is_intron_only(self):
    #     """Checks if variant is entirely intronic"""
    #
    #     # case 1: same base,same sign on offset -> start, end are anchored from the same base
    #     same_start = self._var_c.posedit.pos.start.base == self._var_c.posedit.pos.end.base
    #     pos_offset = self._var_c.posedit.pos.start.offset > 0 and self._var_c.posedit.pos.end.offset > 0
    #     neg_offset = self._var_c.posedit.pos.start.offset < 0 and self._var_c.posedit.pos.end.offset < 0
    #
    #     # case 2: start, end between 2 different exon bases but don't overlap either
    #     same_start_plus_1 = (self._var_c.posedit.pos.start.base + 1) == self._var_c.posedit.pos.end.base
    #     pm_offset = self._var_c.posedit.pos.start.offset > 0 and self._var_c.posedit.pos.end.offset < 0
    #
    #     # other types of introns would overlap at least 1 exon base
    #
    #     return (same_start and (pos_offset or neg_offset)) or (same_start_plus_1 and pm_offset)

    def _incorporate_delins(self):
        """Incorporate delins"""
        seq, cds_start, cds_stop, start, end = self._setup_incorporate()

        ref = self._var_c.posedit.edit.ref
        alt = self._var_c.posedit.edit.alt
        ref_length = end - start if ref is not None else 0    # can't just get from ref since ref isn't always known
        alt_length = len(
            self._var_c.posedit.edit.alt) if self._var_c.posedit.edit.alt is not None else 0
        net_base_change = alt_length - ref_length
        cds_stop += net_base_change

        # incorporate the variant into the sequence (depending on the type)
        is_substitution = False
        if ref is not None and alt is not None:    # delins or SNP
            seq[start:end] = list(alt)
            if len(ref) == 1 and len(alt) == 1:
                is_substitution = True
        elif ref is not None:    # deletion
            del seq[start:end]
        else:    # insertion
            seq[start + 1:start + 1] = list(alt)    # insertion in list before python list index

        if DBG:
            print("net base change: {}".format(net_base_change))
        is_frameshift = net_base_change % 3 != 0
        # use max of mod 3 value and 1 (in event that indel starts in the 5'utr range)
        variant_start_aa = max(int(math.ceil((self._var_c.posedit.pos.start.base) / 3.0)), 1)

        alt_data = AltTranscriptData(
            seq,
            cds_start,
            cds_stop,
            is_frameshift,
            variant_start_aa,
            self._transcript_data.protein_accession,
            is_substitution=is_substitution,
            is_ambiguous=self._ref_has_multiple_stops)
        return alt_data

    def _incorporate_dup(self):
        """Incorporate dup into sequence"""
        seq, cds_start, cds_stop, start, end = self._setup_incorporate()

        dup_seq = seq[start:end]
        seq[end:end] = dup_seq

        is_frameshift = len(dup_seq) % 3 != 0
        variant_start_aa = int(math.ceil((self._var_c.posedit.pos.end.base + 1) / 3.0))

        alt_data = AltTranscriptData(
            seq,
            cds_start,
            cds_stop,
            is_frameshift,
            variant_start_aa,
            self._transcript_data.protein_accession,
            is_ambiguous=self._ref_has_multiple_stops)
        return alt_data

    def _incorporate_inv(self):
        """Incorporate inv into sequence"""
        seq, cds_start, cds_stop, start, end = self._setup_incorporate()

        seq[start:end] = list(reverse_complement(''.join(seq[start:end])))

        is_frameshift = False
        variant_start_aa = max(int(math.ceil((self._var_c.posedit.pos.start.base) / 3.0)), 1)

        alt_data = AltTranscriptData(
            seq,
            cds_start,
            cds_stop,
            is_frameshift,
            variant_start_aa,
            self._transcript_data.protein_accession,
            is_ambiguous=self._ref_has_multiple_stops)
        return alt_data

    def _incorporate_repeat(self):
        """Incorporate repeat int sequence"""
        raise NotImplementedError("hgvs c to p conversion does not support {} type: repeats".format(
            self._var_c))

    def _setup_incorporate(self):
        """Helper to setup incorporate functions
        :return (transcript sequence, cds start [1-based], cds stop [1-based],
        cds start index in seq [inc, 0-based], cds end index in seq [excl, 0-based])
        :rtype (list, int, int, int, int)
        """
        seq = list(self._transcript_data.transcript_sequence)

        # get initial start/end points; will modify these based on the variant length
        cds_start = self._transcript_data.cds_start
        cds_stop = self._transcript_data.cds_stop

        start_end = []
        for pos in (self._var_c.posedit.pos.start, self._var_c.posedit.pos.end):
            # list is zero-based; seq pos is 1-based
            if pos.datum == Datum.CDS_START:
                if pos.base < 0:    # 5' UTR
                    result = cds_start - 1
                else:    # cds/intron
                    if pos.offset <= 0:
                        result = (cds_start - 1) + pos.base - 1
                    else:
                        result = (cds_start - 1) + pos.base
            elif pos.datum == Datum.CDS_END:    # 3' UTR
                result = cds_stop + pos.base - 1
            else:
                raise NotImplementedError("Unsupported/unexpected location")
            start_end.append(result)

        # unpack; increment end by 1 (0-based exclusive)
        (start, end) = start_end
        end += 1

        if DBG:
            print("len seq:{} cds_start:{} cds_stop:{} start:{} end:{}".format(
                len(seq), cds_start, cds_stop, start, end))
        return seq, cds_start, cds_stop, start, end

    def _create_alt_equals_ref_noncds(self):
        """Create an alt seq that matches the reference (for non-cds variants)"""
        alt_data = AltTranscriptData(
            list(self._transcript_data.transcript_sequence),
            self._transcript_data.cds_start,
            self._transcript_data.cds_stop,
            False,
            None,
            self._transcript_data.protein_accession,
            is_ambiguous=True)
        return alt_data

    def _create_no_protein(self):
        """Create a no-protein result"""
        alt_data = AltTranscriptData([],
                                     None,
                                     None,
                                     False,
                                     None,
                                     self._transcript_data.protein_accession,
                                     is_ambiguous=False)
        return alt_data

    def _get_frameshift_start(self, variant_data):
        """Get starting position (AA ref index) of the last frameshift
        which affects the rest of the sequence, i.e. not offset by subsequent frameshifts
        :param variant_data: info on each variant
        :type variant_data: attrs
        :return variant data with additional field for AA index (1-based) of the frameshift start
        :rtype attrs
        """

        if DBG:
            print("is_frameshift:{}".format(variant_data.is_frameshift))
            print("variant_start_aa:{}".format(variant_data.variant_start_aa))
        if variant_data.is_frameshift:
            variant_data.frameshift_start = variant_data.variant_start_aa
        return variant_data


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
