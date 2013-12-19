#
# Tests for conversion of hgvs tags
#
import os
import unittest

import hgvs.hgvsmapper as hgvsmapper
import hgvs.parser
import framework.mock_input_source as mock_input_data_source


class TestHgvsCToP(unittest.TestCase):

    fn = os.path.join(os.path.dirname(__file__), 'data', 'sanity_c_to_p.tsv')
    _datasource = mock_input_data_source.MockInputSource(fn)
    _mapper = hgvsmapper.HGVSMapper(_datasource)
    _parser = hgvs.parser.Parser()

    def test_silent(self):
        hgvsc = "NM_999999.1:c.6A>G"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.="
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_substitution(self):
        hgvsc = "NM_999999.1:c.6A>T"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Lys2Asn"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_substitution_introduces_stop_codon(self):
        hgvsc = "NM_999996.1:c.8C>A"
        hgvsp_expected = "SEGUID_KsJeNmGxgH1t8Wm1Wawsw971Sng:p.Ser3Ter"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_substitution_removes_stop_codon(self):
        hgvsc = "NM_999998.1:c.30G>T"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Ter10Tyrext*3"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_substitution_removes_start_codon(self):
        hgvsc = "NM_999999.1:c.1A>G"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Met1?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_insertion_no_frameshift(self):
        hgvsc = "NM_999999.1:c.6_7insGGG"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Lys2_Ala3insGly"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_insertion_frameshift(self):
        hgvsc = "NM_999999.1:c.22_23insT"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Ala8Valfs*?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_insertion_adds_stop(self):
        hgvsc = "NM_999999.1:c.8_9insTT"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Lys4Ter"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion_no_frameshift(self):
        hgvsc = "NM_999999.1:c.10_12del"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Lys4del"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion2_no_frameshift(self):
        hgvsc = "NM_999999.1:c.4_15del"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Lys2_Ala5del"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion3_no_frameshift_c_term(self):
        hgvsc = "NM_999995.1:c.4_6del"
        hgvsp_expected = "SEGUID_xtu8rRkrlMg5V8HjwWLJFTvJInM:p.Lys3del"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion4_no_frameshift_c_term(self):
        hgvsc = "NM_999994.1:c.4_9del"
        hgvsp_expected = "SEGUID_X0heK3jhOwDrvKEBFMlHupn+lc4:p.Lys3_Lys4del"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion5_no_frameshift(self):
        hgvsc = "NM_999994.1:c.20_25del"
        hgvsp_expected = "SEGUID_X0heK3jhOwDrvKEBFMlHupn+lc4:p.Ala7_Arg9delinsGly"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion6_no_frameshift(self):
        hgvsc = "NM_999999.1:c.5_7del"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Lys2_Ala3delinsThr"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion7_no_frameshift(self):
        hgvsc = "NM_999993.1:c.13_24del"
        hgvsp_expected = "SEGUID_b4UTGUQBBXsyHbKd7Jo+G/HBEh8:p.Arg5_Ala8del"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion_frameshift_nostop(self):
        hgvsc = "NM_999999.1:c.11_12del"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Lys4Serfs*?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion_frameshift_adds_stop(self):
        hgvsc = "NM_999997.1:c.7del"
        hgvsp_expected = "SEGUID_VdwEf25bsURyG6rUUybuNfhSYUc:p.Ala3Argfs*6"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_indel_no_frameshift(self):
        hgvsc = "NM_999999.1:c.11_12delinsTCCCA"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Lys4delinsIlePro"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_indel2_no_frameshift(self):
        hgvsc = "NM_999999.1:c.11_18delinsTCCCA"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Lys4_Phe6delinsIlePro"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_indel_frameshift_nostop(self):
        hgvsc = "NM_999999.1:c.8delinsGG"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Ala3Glyfs*?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_dup_1AA_no_frameshift_2(self):
        hgvsc = "NM_999999.1:c.10_12dup"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Lys4dup"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_dup_1AA_no_frameshift(self):
        hgvsc = "NM_999999.1:c.16_18dup"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Phe6dup"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_dup_2AA_no_frameshift(self):
        hgvsc = "NM_999999.1:c.16_21dup"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Phe6_Arg7dup"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_dup_2AA2_no_frameshift(self):
        hgvsc = "NM_999995.1:c.4_6dup"
        hgvsp_expected = "SEGUID_xtu8rRkrlMg5V8HjwWLJFTvJInM:p.Lys3dup"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_dup_3AA_no_frameshift(self):
        hgvsc = "NM_999999.1:c.16_24dup"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Phe6_Ala8dup"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_dup_frameshift(self):
        hgvsc = "NM_999999.1:c.12_13dup"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.Ala5Glufs*?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_intron(self):
        hgvsc = "NM_999999.1:c.12+1G>A"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_five_prime_utr(self):
        hgvsc = "NM_999999.1:c.-2A>G"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_three_prime_utr(self):
        hgvsc = "NM_999999.1:c.*3G>A"
        hgvsp_expected = "SEGUID_5gk7JyE2f0rvKIE7ZhEDDLBpN9s:p.?"
        self._run_conversion(hgvsc, hgvsp_expected)

    # The following are unsupported
    #
    # def test_repeats(self):
    #     hgvsc = "NM_999999.1:c.12_13[3]"
    #     hgvsp_expected = ""
    #     self._run_conversion(hgvsc, hgvsp_expected)
    #
    # def test_variable_repeats(self):
    #     pass
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

    def _run_conversion(self, hgvsc, hgvsp_expected):
        """Helper method to actually run the test
        :param hgvsc tag
        """
        var_c = TestHgvsCToP._parser.parse_hgvs_variant(hgvsc)
        hgvsp_actual = str(TestHgvsCToP._mapper.hgvsc_to_hgvsp(var_c))
        msg = "hgvsc: {} hgvsp expected: {} actual: {}".format(hgvsc, hgvsp_expected, hgvsp_actual)
        self.assertEqual(hgvsp_expected, hgvsp_actual, msg)

    # TODO - review other classes of hgvs tags (e.g. utr, intronic) - more use cases?
    # 5'utr
    # intronic
    # after stop codon
    # uncertainties in dups/dels (i.e. hgvs tags with ?)



if __name__ == '__main__':
    unittest.main()