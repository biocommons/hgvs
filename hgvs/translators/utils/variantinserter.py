#
# Utility to insert an hgvs variant into a DNA sequence
#
import math

from Bio.Seq import Seq

import hgvs.edit

DBG = True

# TODO - convert the inputs/outputs to whatever the canonical data representation of sequences is

class VariantInserter(object):

    def __init__(self, variant, transcript_data):
        """Constructor

        :param variant representation of hgvs variant
        :type variant
        :param transcript_data representation of transcript
        :type dictionary

        """
        self._variant = variant
        self._transcript_data = transcript_data


    def insert_variant(self):
        """given a variant and a sequence, incorporate the variant and return the new sequence

        Data structure returned is analogous to the data structure used to return the variant sequence,
        but with an additional parameter denoting the start of a frameshift that should affect all bases
        downstream.

        :returns variant sequence data
        :type list of dictionaries
        """

        type_map = {hgvs.edit.RefAlt: self._incorporate_delins,
                    hgvs.edit.Dup: self._incorporate_dup,
                    hgvs.edit.Repeat: self._incorporate_repeat}


        # TODO - loop over each allele rather than assume only 1 variant; return a list for now
        variant_data = []

        has_5utr_variant = False
        has_3utr_variant = False
        has_intron_variant = False

        variant_location = self._get_variant_location()

        # TODO - implement handling non-exons
        if variant_location == "exon":
            edit_type = type(self._variant.posedit.edit)
            this_variant_data = type_map[edit_type]()
        elif variant_location == "intron":
            has_intron_variant = True
        elif variant_location == "three_utr":
            has_3utr_variant = True
        elif variant_location == "five_utr":
            has_5utr_variant = True
        else:   # should never get here
            raise ValueError("value_location = {}".format(variant_location))


        # get the start of the "terminal" frameshift (i.e. one never "cancelled out")
        this_variant_data = self._get_frameshift_start(this_variant_data)

        variant_data.append(this_variant_data)

        return variant_data


    def _get_variant_location(self):
        """Categorize variant by location in transcript (5'utr, exon, intron, 3'utr)
        :return "exon", "intron", "five_utr", "three_utr"
        """
        # TODO - implement
        return "exon"


    def _incorporate_delins(self):
        """Incorporate delins
        """
        seq = list(self._transcript_data['transcript_sequence'])

        # get initial start/end points; will modify these based on the variant length
        cds_start = self._transcript_data['cds_start']
        cds_stop = self._transcript_data['cds_stop']

        start = (cds_start - 1) + self._variant.posedit.pos.start.base - 1   # list is zero-based; seq pos is 1-based
        end = (cds_start - 1) + self._variant.posedit.pos.end.base
        ref = self._variant.posedit.edit.ref
        alt = self._variant.posedit.edit.alt
        ref_length = end - start

        if DBG:
            print "start: {} end:{} ref:{} alt{}".format(start, end, ref, alt)

        if ref is not None and alt is not None: # delins (or snp)
            seq[start:end] = list(alt)
            if len(alt) != len(ref):
                cds_stop += (len(alt) - (end - start))
            net_base_change = len(alt) - ref_length
        elif ref is not None:       # deletion
            del seq[start:end]
            cds_stop -= end - start
            net_base_change = -ref_length
        else:           # insertion (only alt)
            seq[start + 1:start + 1] = list(alt)    # insertion in list before python list index
            cds_stop += len(alt)
            net_base_change = len(alt)

        if DBG:
            print "Net base change: {}".format(net_base_change)
        is_frameshift = net_base_change % 3 != 0
        AA_start = int(math.ceil((self._variant.posedit.pos.start.base + 1) / 3.0))

        variant_data = self._create_AA_variant_output(seq, cds_start, cds_stop, is_frameshift, AA_start,
                                                      self._transcript_data['protein_accession'])
        return variant_data


    def _incorporate_dup(self):
        """Incorporate dup into sequence
        """
        seq = list(self._transcript_data['transcript_sequence'])

        # get initial start/end points; will modify these based on the variant length
        cds_start = self._transcript_data['cds_start']
        cds_stop = self._transcript_data['cds_stop']

        start = (cds_start - 1) + self._variant.posedit.pos.start.base - 1   # list is zero-based; seq pos is 1-based
        end = (cds_start - 1) + self._variant.posedit.pos.end.base

        dup_seq = seq[start:end]

        seq[end:end] = dup_seq

        is_frameshift = len(dup_seq) % 3 != 0
        AA_start = int(math.ceil((self._variant.posedit.pos.end.base + 1) / 3.0))

        variant_data = self._create_AA_variant_output(seq, cds_start, cds_stop, is_frameshift, AA_start,
                                                      self._transcript_data['protein_accession'])
        return variant_data

    def _incorporate_repeat(self):
        """Incorporate repeat int sequence
        """
        # TODO - implement
        return self._transcript_data


    def _get_frameshift_start(self, variant_data):
        """Get starting position (AA ref index) of the last frameshift
        which affects the rest of the sequence, i.e. not offset by subsequent frameshifts
        :param variant_data: info on each variant
        :type list of dictionaries
        :return variant data with additional field for AA index (1-based) of the frameshift start; -1 if none
        """
        # TODO - implement for 2+ variants

        if variant_data['is_frameshift']:
            frameshift_start = variant_data['AA_start']
        else:
            frameshift_start = -1
        variant_data['frameshift_start'] = frameshift_start
        return variant_data


    def _create_AA_variant_output(self, seq, cds_start, cds_stop, is_frameshift, AA_start, accession):
        """Common code for creating a variant sequence in the proper format once the sequence has been modified
        :param seq: DNA sequence wiith variant incorporated
        :type str
        :param cds_start: coding sequence start (1-based)
        :type int
        :param cds_stop: coding sequence stop (1-based)
        :type int
        :param is_frameshift: is this variant a frameshift
        :type bool
        :param AA_start: AA start index (1-based) for this variant
        :type int
        :param accession: protein accession, e.g. NP_999999.2
        :type str
        :return variant sequence data
        :type dict
        """

        # padding list so biopython won't complain during the conversion
        if len(seq) % 3 != 0:
            seq.extend(['N']*(3-len(seq) % 3))

        seq = ''.join(seq)

        seq_cds = Seq(seq[cds_start - 1:])
        seq_AA = str(seq_cds.translate())
        stop_pos = seq_AA.find("*")
        if stop_pos != -1:
            seq_AA = seq_AA[:stop_pos + 1]

        variant_data = {
            'transcript_sequence': seq,
            'AA_sequence': seq_AA,
            'cds_start': cds_start,
            'cds_stop': cds_stop,
            'is_frameshift': is_frameshift,
            'AA_start': AA_start,
            'protein_accession': accession
        }

        return variant_data

def main():
    pass


if __name__ == "__main__":
    main()

