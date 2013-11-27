#
# Tests for conversion of hgvs tags
#

import unittest

import hgvs.translators.hgvs_c_to_p as hgvs_c_to_p
import framework.mock_input_source as mock_input_data_source

class TestHgvsCToP(unittest.TestCase):


    _data_source = mock_input_data_source.MockInputSource('data/transcript_data.tsv')
    _translator = hgvs_c_to_p.HgvsCToP(_data_source)

    def test_silent(self):
        hgvsc = "NM_999999.1:c.6A>G"
        expected_hgvsp = "NP_999999.1:p.="
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_substitution(self):
        hgvsc = "NM_999999.1:c.6A>T"
        expected_hgvsp = "NP_999999.1:p.Lys2Asn"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_substitution_introduces_stop_codon(self):
        hgvsc = "NM_999996.1:c.8C>A"
        expected_hgvsp = "NP_999996.1:p.Ser3Ter"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_substitution_removes_stop_codon(self):
        hgvsc = "NM_999998.1:c.30G>T"
        expected_hgvsp = "NP_999998.1:p.Ter10Argext*3"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_insertion_no_frameshift(self):
        hgvsc = "NM_999999.1:c.6_7insGGG"
        expected_hgvsp = "NP_999999.1:p.Ala3_Lys4insGly"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_insertion_frameshift(self):
        hgvsc = "NM_999999.1:c.22_23insT"
        expected_hgvsp = "NP_999999.1:p.Ala8Valfs*?"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_insertion_adds_stop(self):
        hgvsc = "NM_999999.1:c.8_9insTT"
        expected_hgvsp = "NP_999999.1:p.Lys4Ter"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_deletion_no_frameshift(self):
        hgvsc = "NM_999999.1:c.10_12del"
        expected_hgvsp = "NP_999999.1:p.Lys4del"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_deletion2_no_frameshift(self):
        hgvsc = "NM_999999.1:c.4_15del"
        expected_hgvsp = "NP_999999.1:p.Lys2_Ala5del"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_deletion_frameshift_nostop(self):
        hgvsc = "NM_999999.1:c.11_12del"
        expected_hgvsp = "NP_999999.1:p.Lys4Serfs*?"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_deletion_frameshift_adds_stop(self):
        hgvsc = "NM_999997.1:c.7del"
        expected_hgvsp = "NP_999997.1:p.Ala3Argfs*6"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_indel_no_frameshift(self):
        hgvsc = "NM_999999.1:c.11_12delinsTCCCA"
        expected_hgvsp = "NP_999999.1:p.Lys4delinsIlePro"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_indel2_no_frameshift(self):
        hgvsc = "NM_999999.1:c.11_18delinsTCCCA"
        expected_hgvsp = "NP_999999.1:p.Lys4_Phe6delinsIlePro"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_indel_frameshift_nostop(self):
        hgvsc = "NM_999999.1:c.8delinsGG"
        expected_hgvsp = "NP_999999.1:p.Ala3Glyfs*?"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_dup_1AA_no_frameshift_2(self):
        hgvsc = "NM_999999.1:c.10_12dup"
        expected_hgvsp = "NP_999999.1:p.Lys4dup"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_dup_1AA_no_frameshift(self):
        hgvsc = "NM_999999.1:c.16_18dup"
        expected_hgvsp = "NP_999999.1:p.Phe6dup"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_dup_2AA_no_frameshift(self):
        hgvsc = "NM_999999.1:c.16_21dup"
        expected_hgvsp = "NP_999999.1:p.Phe6_Arg7dup"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_dup_3AA_no_frameshift(self):
        hgvsc = "NM_999999.1:c.16_24dup"
        expected_hgvsp = "NP_999999.1:p.Phe6_Ala8dup"
        self._run_conversion(hgvsc, expected_hgvsp)

    def test_dup_frameshift(self):
        hgvsc = "NM_999999.1:c.12_13dup"
        expected_hgvsp = "NP_999999.1:p.Ala5Glufs*?"
        self._run_conversion(hgvsc, expected_hgvsp)

    # The following are unsupported
    #
    # def test_repeats(self):
    #     hgvsc = "NM_999999.1:c.12_13[3]"
    #     expected_hgvsp = ""
    #     self._run_conversion(hgvsc, expected_hgvsp)
    #
    # def test_variable_repeats(self):
    #     pass
    #
    # def test_extension_new_translation_initiation_site(self):
    #     hgvsc = "NM_999999.1:c.1A>G"
    #     expected_hgvsp = ""
    #     self._run_conversion(hgvsc, expected_hgvsp)
    #
    # def test_indeterminate_entire_exon_del(self):
    #     pass
    #
    # def test_indeterminate_entire_exon_dup(self):
    #     pass
    #
    # def test_mosaic(self):
    #     pass
    #
    # def test_chimera(self):
    #     pass
    #
    # def test_two_changes_same_allele(self):
    #     pass
    #
    # def test_two_changes_diff_allele(self):
    #     pass
    #
    # def test_two_changes_unknown_allele(self):
    #     pass

    def _run_conversion(self, hgvsc, expected_hgvsp):
        """Helper method to actually run the test
        :param hgvsc tag
        """
        actual_hgvsp = str(TestHgvsCToP._translator.convert(hgvsc))
        msg = "hgvsp expected: {} actual: {}".format(expected_hgvsp, actual_hgvsp)
        self.assertEqual(expected_hgvsp, actual_hgvsp, msg)

    # TODO - review other classes of hgvs tags (e.g. utr, intronic) - more use cases?
    # 5'utr
    # intronic
    # after stop codon
    # uncertainties in dups/dels (i.e. hgvs tags with ?)



if __name__ == '__main__':
    unittest.main()
