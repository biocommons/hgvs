import unittest

import hgvs.parser
import hgvs.translators.utils.variantinserter as variantinserter

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

    # def test_2_substitutions(self):
    #     pass
    #
    # def test_2_indel_no_net_frameshift(self):
    #     pass
    #
    # def test_2_indel_net_frameshift(self):
    #     pass

    def _run_comparison(self, hgvsc, expected_sequence):
        var = self._parser.hgvs_variant(hgvsc)
        transcript_data = self._datasource.get_sequence(var.seqref)
        inserter = variantinserter.VariantInserter(var, transcript_data)
        insert_result = inserter.insert_variant()
        actual_sequence = insert_result[0]['transcript_sequence']
        msg = "expected: {}\nactual  : {}".format(expected_sequence, actual_sequence)
        self.assertEqual(expected_sequence, actual_sequence, msg)


if __name__ == '__main__':
    unittest.main()
