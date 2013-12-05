import unittest

from Bio.Seq import Seq

import hgvs.parser
import hgvs.translators.hgvs_c_to_p as hgvs_c_to_p
import hgvs.translators.tools.altseqbuilder as altseqbuilder

import framework.mock_input_source as mock_input_data_source

class TestVariantInserter(unittest.TestCase):


    # root sequence = ""
    _datasource = mock_input_data_source.MockInputSource('data/transcript_data.tsv')
    _parser = hgvs.parser.Parser()


    def test_substitution_start(self):
        hgvsc = "NM_999999.1:c.1A>T"
        expected_sequence = "AAAATCAAATTGAAAGCGAAAGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_substitution_middle(self):
        hgvsc = "NM_999999.1:c.6A>T"
        expected_sequence = "AAAATCAAAATGAATGCGAAAGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_substitution_end(self):
        hgvsc = "NM_999999.1:c.30G>C"
        expected_sequence = "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGCGAAATACGGG"
        self._run_comparison(hgvsc, expected_sequence)

    # TODO - build in support when system can handle variants in 5'utr region
    # def test_insertion_before_start(self):
    #     hgvsc = "NM_999999.1:c.-1_1insGGG"
    #     expected_sequence = "AAAATCAAAGGGATGAAAGCGAAAGCGTTTCGCGCGAAATAGGGG"
    #     self._run_comparison(hgvsc, expected_sequence)

    def test_insertion_start(self):
        hgvsc = "NM_999999.1:c.1_2insAAA"
        expected_sequence = "AAAATCAAAAAAATGAAAGCGAAAGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_insertion_middle(self):
        hgvsc = "NM_999999.1:c.22_23insT"
        expected_sequence = "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGTCGAAATAGGGGNN"
        self._run_comparison(hgvsc, expected_sequence)

    def test_insertion_end(self):
        hgvsc = "NM_999999.1:c.29_30insGG"
        expected_sequence = "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGCGAAATAGGGGGGN"
        self._run_comparison(hgvsc, expected_sequence)

    # TODO - build in support when system can handle variants in 3'utr region
    # def test_insertion_after_end(self):
    #     hgvsc = "NM_999999.1:c.30_*1insAA"
    #     expected_sequence = "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGCGAAATAGAAGGGN"
    #     self._run_comparison(hgvsc, expected_sequence)

    def test_deletion_start(self):
        hgvsc = "NM_999999.1:c.1del"
        expected_sequence = "AAAATCAAATGAAAGCGAAAGCGTTTCGCGCGAAATAGGGGN"
        self._run_comparison(hgvsc, expected_sequence)

    def test_deletion_middle(self):
        hgvsc = "NM_999999.1:c.2_7del"
        expected_sequence = "AAAATCAAAACGAAAGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_deletion_end(self):
        hgvsc = "NM_999999.1:c.30del"
        expected_sequence = "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGCGAAATAGGGN"
        self._run_comparison(hgvsc, expected_sequence)

    def test_delins_start(self):
        hgvsc = "NM_999999.1:c.1delinsTTTT"
        expected_sequence = "AAAATCAAATTTTTGAAAGCGAAAGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_delins_middle(self):
        hgvsc = "NM_999999.1:c.2_3delinsAA"
        expected_sequence = "AAAATCAAAAAAAAAGCGAAAGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_delins_end(self):
        hgvsc = "NM_999999.1:c.30delinsCCCC"
        expected_sequence = "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGCGAAATACCCCGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_dup(self):
        hgvsc = "NM_999999.1:c.16_24dup"
        expected_sequence = "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)


    # def test_2_substitutions(self):
    #     pass
    #
    # def test_2_indel_no_net_frameshift(self):
    #     pass
    #
    # def test_2_indel_net_frameshift(self):
    #     pass

    def _run_comparison(self, hgvsc, expected_sequence):


        var = self._parser.parse_hgvs_variant(hgvsc)

        ts_exons = self._datasource.fetch_transcript_exons(var.seqref, '_')
        ts_info = self._datasource.fetch_transcript_info(var.seqref)

        cds_start = ts_info['cds_start_i'] + 1
        cds_stop = ts_info['cds_stop_i']

        # concatenate sequences by exon
        ord_seq = {x['ord']: x['t_seq_a'] for x in ts_exons}
        ordered_keys = ord_seq.keys()
        ordered_keys.sort()
        transcript_sequence = ''.join(ord_seq[x] for x in ordered_keys)

        transcript_dna_cds = Seq(transcript_sequence[cds_start - 1:cds_stop])
        protein_seq = str(transcript_dna_cds.translate())
        protein_acc = hgvs.stopgap.pseq_to_ac(protein_seq)

        transcript_data = hgvs_c_to_p.RefTranscriptData(transcript_sequence,
                                            protein_seq,
                                            cds_start,
                                            cds_stop,
                                            protein_acc)

        builder = altseqbuilder.AltSeqBuilder(var, transcript_data)
        insert_result = builder.build_altseq()
        actual_sequence = insert_result[0].transcript_sequence
        msg = "expected: {}\nactual  : {}".format(expected_sequence, actual_sequence)
        self.assertEqual(expected_sequence, actual_sequence, msg)


if __name__ == '__main__':
    unittest.main()
