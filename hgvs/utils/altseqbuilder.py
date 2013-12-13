#
# Utility to insert an hgvs variant into a transcript sequence.
# Generates a record corresponding to the modified transcript sequence,
# along with annotations for use in conversion to an hgvsp tag.
# Used in hgvsc to hgvsp conversion.
#
import math
import recordtype

from Bio.Seq import Seq

import hgvs.edit
from hgvs.location import CDS_END

DBG = False

class AltTranscriptData(recordtype.recordtype('AltTranscriptData', [
        'transcript_sequence', 'aa_sequence', 'cds_start', 'cds_stop', 'protein_accession',
        ('is_frameshift', False), ('variant_start_aa', None), ('frameshift_start', None), ('is_substitution', False)])):

    @classmethod
    def create_for_variant_inserter(cls, seq, cds_start, cds_stop, is_frameshift, variant_start_aa, accession,
                                    pad_seq = True, is_substitution=False):
        """Create a variant sequence using inputs from VariantInserter
        :param seq: DNA sequence wiith variant incorporated
        :type str
        :param cds_start: coding sequence start (1-based)
        :type int
        :param cds_stop: coding sequence stop (1-based)
        :type int
        :param is_frameshift: is this variant a frameshift
        :type bool
        :param variant_start_aa: AA start index (1-based) for this variant
        :type int
        :param accession: protein accession, e.g. NP_999999.2
        :type str
        :param pad_seq: if true, pad sequence to be divisible by three (default true)
        :type bool
        :return variant sequence data
        :type recordtype
        """

        # padding list so biopython won't complain during the conversion
        if pad_seq and len(seq) % 3 != 0:
            seq.extend(['N']*(3-len(seq) % 3))

        seq = ''.join(seq)

        seq_cds = Seq(seq[cds_start - 1:])
        seq_aa = str(seq_cds.translate())
        stop_pos = seq_aa.find("*")
        if stop_pos != -1:
            seq_aa = seq_aa[:stop_pos + 1]

        alt_data = AltTranscriptData(seq, seq_aa, cds_start, cds_stop, accession, is_frameshift, variant_start_aa,
                                     is_substitution=is_substitution)

        return alt_data


class AltSeqBuilder(object):

    EXON = "exon"
    INTRON = "intron"
    F_UTR = "five utr"
    T_UTR = "three utr"

    def __init__(self, var_c, transcript_data):
        """Constructor

        :param variant representation of hgvs variant
        :type variant
        :param transcript_data representation of transcript
        :type recordtype

        """
        self._var_c = var_c
        self._transcript_data = transcript_data
        if DBG:
            print transcript_data.transcript_sequence

    def build_altseq(self):
        """given a variant and a sequence, incorporate the variant and return the new sequence

        Data structure returned is analogous to the data structure used to return the variant sequence,
        but with an additional parameter denoting the start of a frameshift that should affect all bases
        downstream.

        :returns variant sequence data
        :type list of dictionaries
        """
        NOT_CDS = "not_cds_variant"

        type_map = {hgvs.edit.NARefAlt: self._incorporate_delins,
                    hgvs.edit.Dup: self._incorporate_dup,
                    hgvs.edit.Repeat: self._incorporate_repeat,
                    NOT_CDS: self._create_alt_eq_ref
                    }

        # TODO - loop over each allele rather than assume only 1 variant; return a list for now
        alt_data = []

        variant_location = self._get_variant_location()

        if variant_location == self.EXON:
            edit_type = type(self._var_c.posedit.edit)
        elif variant_location == self.INTRON:
            edit_type = NOT_CDS
        elif variant_location == self.T_UTR:
            edit_type = NOT_CDS
        elif variant_location == self.F_UTR:
            # TODO - handle case where variant introduces a Met (new start)
            edit_type = NOT_CDS
        else:   # should never get here
            raise ValueError("value_location = {}".format(variant_location))

        try:
            this_alt_data = type_map[edit_type]()
        except KeyError as e:
            raise NotImplementedError("c to p translation unsupported for {} type {}".format(self._var_c, edit_type))

        # get the start of the "terminal" frameshift (i.e. one never "cancelled out")
        this_alt_data = self._get_frameshift_start(this_alt_data)
        alt_data.append(this_alt_data)
        if DBG:
            print this_alt_data.transcript_sequence

        return alt_data

    def _get_variant_location(self):
        """Categorize variant by location in transcript (5'utr, exon, intron, 3'utr)
        :return "exon", "intron", "five_utr", "three_utr"
        """
        if CDS_END in [self._var_c.posedit.pos.start.datum, self._var_c.posedit.pos.end.datum]:
            result = self.T_UTR
        elif self._var_c.posedit.pos.start.base < 0 or self._var_c.posedit.pos.end.base < 0:
            result = self.F_UTR
        elif self._var_c.posedit.pos.start.offset != 0 or self._var_c.posedit.pos.end.offset != 0:
            result = self.INTRON
        else:
            result = self.EXON
        return result

    def _incorporate_delins(self):
        """Incorporate delins"""
        seq, cds_start, cds_stop, start, end = self._setup_incorporate()

        ref = self._var_c.posedit.edit.ref
        alt = self._var_c.posedit.edit.alt
        ref_length = end - start if ref is not None else 0  # can't just get from ref since ref isn't always known
        alt_length = len(self._var_c.posedit.edit.alt) if self._var_c.posedit.edit.alt is not None else 0
        net_base_change = alt_length - ref_length
        cds_stop += net_base_change

        # incorporate the variant into the sequence (depending on the type)
        is_substitution = False
        if ref is not None and alt is not None:     # delins or SNP
            seq[start:end] = list(alt)
            if len(ref) == 1 and len(alt) == 1:
                is_substitution = True
        elif ref is not None:                       # deletion
            del seq[start:end]
        else:                                       # insertion
            seq[start + 1:start + 1] = list(alt)    # insertion in list before python list index

        is_frameshift = net_base_change % 3 != 0
        variant_start_aa = int(math.ceil((self._var_c.posedit.pos.start.base) / 3.0))

        alt_data = AltTranscriptData.create_for_variant_inserter(seq, cds_start, cds_stop,
                                                               is_frameshift, variant_start_aa,
                                                               self._transcript_data.protein_accession,
                                                               is_substitution=is_substitution)
        return alt_data

    def _incorporate_dup(self):
        """Incorporate dup into sequence"""
        seq, cds_start, cds_stop, start, end = self._setup_incorporate()

        dup_seq = seq[start:end]
        seq[end:end] = dup_seq
        
        is_frameshift = len(dup_seq) % 3 != 0
        variant_start_aa = int(math.ceil((self._var_c.posedit.pos.end.base + 1) / 3.0))

        alt_data = AltTranscriptData.create_for_variant_inserter(seq, cds_start, cds_stop,
                                                                 is_frameshift, variant_start_aa,
                                                                 self._transcript_data.protein_accession)
        return alt_data

    def _incorporate_repeat(self):
        """Incorporate repeat int sequence"""
        # TODO - implement
        raise NotImplementedError("hgvs c to p conversion does not support {} type: repeats".format(self._var_c))

    def _setup_incorporate(self):
        """Helper to setup incorporate functions
        :return (transcript sequence, cds start [1-based], cds stop [1-based],
        cds start index in seq [inc, 0-based], cds end index in seq [excl, 0-based])
        :type (str, int, int, int, int)
        """
        seq = list(self._transcript_data.transcript_sequence)

        # get initial start/end points; will modify these based on the variant length
        cds_start = self._transcript_data.cds_start
        cds_stop = self._transcript_data.cds_stop

        start = (cds_start - 1) + self._var_c.posedit.pos.start.base - 1   # list is zero-based; seq pos is 1-based
        end = (cds_start - 1) + self._var_c.posedit.pos.end.base

        return seq, cds_start, cds_stop, start, end

    def _create_alt_eq_ref(self):
        """Create an alt seq that matches the reference"""
        alt_data = AltTranscriptData.create_for_variant_inserter(self._transcript_data.transcript_sequence,
                                                                 self._transcript_data.cds_start,
                                                                 self._transcript_data.cds_stop,
                                                                 False,
                                                                 None,
                                                                 self._transcript_data.protein_accession,
                                                                 pad_seq=False)
        return alt_data

    def _get_frameshift_start(self, variant_data):
        """Get starting position (AA ref index) of the last frameshift
        which affects the rest of the sequence, i.e. not offset by subsequent frameshifts
        :param variant_data: info on each variant
        :type list of dictionaries
        :return variant data with additional field for AA index (1-based) of the frameshift start; -1 if none
        """
        # TODO - implement for 2+ variants

        if variant_data.is_frameshift:
            variant_data.frameshift_start = variant_data.variant_start_aa
        return variant_data