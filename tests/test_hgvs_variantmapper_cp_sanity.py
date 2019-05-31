# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

#
# Tests for conversion of hgvs tags
#
import os
import unittest

import hgvs.variantmapper as variantmapper
import hgvs.parser

import support.mock_input_source as mock_input_data_source


class TestHgvsCToP(unittest.TestCase):

    fn = os.path.join(os.path.dirname(__file__), "data", "sanity_cp.tsv")
    _datasource = mock_input_data_source.MockInputSource(fn)
    _mapper = variantmapper.VariantMapper(_datasource, prevalidation_level="INTRINSIC")
    _parser = hgvs.parser.Parser()

    def test_silent(self):
        hgvsc = "NM_999999.1:c.6A>G"
        hgvsp_expected = "MOCK:p.(Lys2=)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_substitution(self):
        hgvsc = "NM_999999.1:c.6A>T"
        hgvsp_expected = "MOCK:p.(Lys2Asn)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_substitution_introduces_stop_codon(self):
        hgvsc = "NM_999996.1:c.8C>A"
        hgvsp_expected = "MOCK:p.(Ser3Ter)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_substitution_removes_stop_codon(self):
        hgvsc = "NM_999998.1:c.30G>T"
        hgvsp_expected = "MOCK:p.(Ter10TyrextTer3)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_insertion_no_frameshift(self):
        hgvsc = "NM_999999.1:c.6_7insGGG"
        hgvsp_expected = "MOCK:p.(Lys2_Ala3insGly)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_insertion_frameshift(self):
        hgvsc = "NM_999999.1:c.22_23insT"
        hgvsp_expected = "MOCK:p.(Ala8ValfsTer?)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_insertion_adds_stop(self):
        hgvsc = "NM_999999.1:c.8_9insTT"
        hgvsp_expected = "MOCK:p.(Lys4Ter)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion_no_frameshift(self):
        hgvsc = "NM_999999.1:c.10_12del"
        hgvsp_expected = "MOCK:p.(Lys4del)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion2_no_frameshift(self):
        hgvsc = "NM_999999.1:c.4_15del"
        hgvsp_expected = "MOCK:p.(Lys2_Ala5del)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion3_no_frameshift_c_term(self):
        hgvsc = "NM_999995.1:c.4_6del"
        hgvsp_expected = "MOCK:p.(Lys3del)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion4_no_frameshift_c_term(self):
        hgvsc = "NM_999994.1:c.4_9del"
        hgvsp_expected = "MOCK:p.(Lys3_Lys4del)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion5_no_frameshift(self):
        hgvsc = "NM_999994.1:c.20_25del"
        hgvsp_expected = "MOCK:p.(Ala7_Arg9delinsGly)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion6_no_frameshift(self):
        hgvsc = "NM_999999.1:c.5_7del"
        hgvsp_expected = "MOCK:p.(Lys2_Ala3delinsThr)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion7_no_frameshift(self):
        hgvsc = "NM_999993.1:c.13_24del"
        hgvsp_expected = "MOCK:p.(Arg5_Ala8del)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion_frameshift_nostop(self):
        hgvsc = "NM_999999.1:c.11_12del"
        hgvsp_expected = "MOCK:p.(Lys4SerfsTer?)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion_frameshift_adds_stop(self):
        hgvsc = "NM_999997.1:c.7del"
        hgvsp_expected = "MOCK:p.(Ala3ArgfsTer6)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion_no_frameshift_removes_stop_plus_previous(self):
        hgvsc = "NM_999999.1:c.25_30del"
        hgvsp_expected = "MOCK:p.(Lys9_Ter10delinsGly)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_indel_no_frameshift(self):
        hgvsc = "NM_999999.1:c.11_12delinsTCCCA"
        hgvsp_expected = "MOCK:p.(Lys4delinsIlePro)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_indel2_no_frameshift(self):
        hgvsc = "NM_999999.1:c.11_18delinsTCCCA"
        hgvsp_expected = "MOCK:p.(Lys4_Phe6delinsIlePro)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_indel_frameshift_nostop(self):
        hgvsc = "NM_999999.1:c.8delinsGG"
        hgvsp_expected = "MOCK:p.(Ala3GlyfsTer?)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_dup_1AA_no_frameshift_2(self):
        hgvsc = "NM_999999.1:c.10_12dup"
        hgvsp_expected = "MOCK:p.(Lys4dup)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_dup_1AA_no_frameshift(self):
        hgvsc = "NM_999999.1:c.16_18dup"
        hgvsp_expected = "MOCK:p.(Phe6dup)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_dup_2AA_no_frameshift(self):
        hgvsc = "NM_999999.1:c.16_21dup"
        hgvsp_expected = "MOCK:p.(Phe6_Arg7dup)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_dup_2AA2_no_frameshift(self):
        hgvsc = "NM_999995.1:c.4_6dup"
        hgvsp_expected = "MOCK:p.(Lys3dup)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_dup_3AA_no_frameshift(self):
        hgvsc = "NM_999999.1:c.16_24dup"
        hgvsp_expected = "MOCK:p.(Phe6_Ala8dup)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_dup_frameshift(self):
        hgvsc = "NM_999999.1:c.12_13dup"
        hgvsp_expected = "MOCK:p.(Ala5GlufsTer?)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_intron(self):
        hgvsc = "NM_999999.1:c.12+1G>A"
        hgvsp_expected = "MOCK:p.?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_five_prime_utr(self):
        hgvsc = "NM_999999.1:c.-2A>G"
        hgvsp_expected = "MOCK:p.?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_three_prime_utr(self):
        hgvsc = "NM_999999.1:c.*3G>A"
        hgvsp_expected = "MOCK:p.?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion_into_three_prime_utr_frameshift(self):
        hgvsc = "NM_999999.1:c.27_*3del"
        hgvsp_expected = "MOCK:p.(Lys9XaafsTer?)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion_into_three_prime_utr_no_frameshift(self):
        hgvsc = "NM_999995.1:c.28_*3del"
        hgvsp_expected = "MOCK:p.(Lys10_Ter11delinsArgGlnPheArg)"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_delins_into_three_prime_utr_no_frameshift(self):
        hgvsc = "NM_999995.1:c.28_*3delinsGGG"
        hgvsp_expected = "MOCK:p.(Lys10_Ter11delinsGlyArgGlnPheArg)"
        self._run_conversion(hgvsc, hgvsp_expected)

    # See recommendations re p.? (p.Met1?) at:
    # http://varnomen.hgvs.org/recommendations/protein/variant/substitution/
    def test_substitution_removes_start_codon(self):
        hgvsc = "NM_999999.1:c.1A>G"
        hgvsp_expected = "MOCK:p.Met1?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion_from_five_prime_utr_frameshift(self):
        hgvsc = "NM_999999.1:c.-3_1del"
        hgvsp_expected = "MOCK:p.Met1?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_deletion_from_five_prime_utr_no_frameshift(self):
        hgvsc = "NM_999999.1:c.-3_3del"
        hgvsp_expected = "MOCK:p.Met1?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_delins_from_five_prime_utr_no_frameshift(self):
        hgvsc = "NM_999999.1:c.-3_3delinsAAA"
        hgvsp_expected = "MOCK:p.Met1?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_delete_entire_gene(self):
        hgvsc = "NM_999999.1:c.-3_*1del"
        hgvsp_expected = "MOCK:p.0?"
        self._run_conversion(hgvsc, hgvsp_expected)

    def test_multiple_stop_codons(self):
        hgvsc = "NM_999992.1:c.4G>A"
        hgvsp_expected = "MOCK:p.?"
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
        ac_p = "MOCK"
        hgvsp_actual = str(TestHgvsCToP._mapper.c_to_p(var_c, ac_p))
        msg = "hgvsc: {} hgvsp expected: {} actual: {}".format(hgvsc, hgvsp_expected, hgvsp_actual)
        self.assertEqual(hgvsp_expected, hgvsp_actual, msg)

    # TODO - review other classes of hgvs tags (e.g. utr, intronic) - more use cases?
    # 5'utr
    # intronic
    # after stop codon
    # uncertainties in dups/dels (i.e. hgvs tags with ?)


if __name__ == "__main__":
    unittest.main()

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
