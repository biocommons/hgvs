# -*- coding: utf-8 -*-
import os
import unittest

import pytest

from support import CACHE

import hgvs
import hgvs.parser
import hgvs.sequencevariant


def test_gene_formatting(parser):
    v = parser.parse("NM_01234.5(BOGUS):c.65A>C")
    assert str(v) == "NM_01234.5(BOGUS):c.65A>C"
    v.gene = None
    assert str(v) == "NM_01234.5:c.65A>C"


@pytest.mark.quick
@pytest.mark.models
class Test_SequenceVariant(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.hdp = hgvs.dataproviders.uta.connect(
            mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE
        )
        cls.vm = hgvs.variantmapper.VariantMapper(cls.hdp)
        cls.hp = hgvs.parser.Parser()


    def test_SequenceVariant(self):
        var = hgvs.sequencevariant.SequenceVariant(ac="AC", type="B", posedit="1234DE>FG")
        self.assertEqual(str(var), "AC:B.1234DE>FG")

    def test_fill_ref(self):

        # fill reference for sequence variants
        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.31_32del").fill_ref(self.hdp)
        self.assertEqual(var.format({"max_ref_length": None}), "NM_001166478.1:c.31_32delTT")

        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.31_32del2").fill_ref(self.hdp)
        self.assertEqual(var.format({"max_ref_length": None}), "NM_001166478.1:c.31_32delTT")

        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.2_7delinsTTTAGA").fill_ref(self.hdp)
        self.assertEqual(
            var.format({"max_ref_length": None}), "NM_001166478.1:c.2_7delTGAAGAinsTTTAGA"
        )

        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.35_36dup").fill_ref(self.hdp)
        self.assertEqual(var.format({"max_ref_length": None}), "NM_001166478.1:c.35_36dupTC")

        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.18_19insACT").fill_ref(self.hdp)
        self.assertEqual(var.format({"max_ref_length": None}), "NM_001166478.1:c.18_19insACT")

        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.31=").fill_ref(self.hdp)
        self.assertEqual(var.format({"max_ref_length": None}), "NM_001166478.1:c.31T=")

    def test_format(self):

        # Global default settings
        var = self.hp.parse_hgvs_variant("NP_001628.1:p.Gly528Arg")
        self.assertEqual(str(var), "NP_001628.1:p.Gly528Arg")
        self.assertEqual(var.format(), "NP_001628.1:p.Gly528Arg")

        # Change global settings
        hgvs.global_config.formatting.p_3_letter = False
        self.assertEqual(str(var), "NP_001628.1:p.G528R")

        # Custom settings
        hgvs.global_config.formatting.p_3_letter = True
        conf = {"p_3_letter": False}
        self.assertEqual(var.format(conf), "NP_001628.1:p.G528R")

        var = self.hp.parse_hgvs_variant("NP_001628.1:p.Gly528Ter")
        conf = {"p_term_asterisk": True}
        self.assertEqual(var.format(conf), "NP_001628.1:p.Gly528*")
        self.assertEqual(var.format(), "NP_001628.1:p.Gly528Ter")
        conf = {"p_3_letter": False}
        self.assertEqual(var.format(conf), "NP_001628.1:p.G528*")
        self.assertEqual(var.format(), "NP_001628.1:p.Gly528Ter")

        # Remove reference sequence
        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.31_32delTT")
        self.assertEqual(str(var), "NM_001166478.1:c.31_32del")
        self.assertEqual(var.format(conf={"max_ref_length": 1}), "NM_001166478.1:c.31_32del")
        self.assertEqual(var.format(conf={"max_ref_length": 2}), "NM_001166478.1:c.31_32delTT")
        self.assertEqual(var.format(conf={"max_ref_length": None}), "NM_001166478.1:c.31_32delTT")

        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.31_32del2")
        self.assertEqual(str(var), "NM_001166478.1:c.31_32del")
        self.assertEqual(var.format(conf={"max_ref_length": None}), "NM_001166478.1:c.31_32del2")

        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.31_32delTTinsAA")
        self.assertEqual(str(var), "NM_001166478.1:c.31_32delinsAA")
        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.35_36dupTC")
        self.assertEqual(str(var), "NM_001166478.1:c.35_36dup")
        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.31T=")
        self.assertEqual(str(var), "NM_001166478.1:c.31=")
        self.assertEqual(var.format(conf={"max_ref_length": None}), "NM_001166478.1:c.31T=")

    def test_uncertain(self):

        vs = "NC_000005.9:g.(90136803_90144453)_(90159675_90261231)dup"
        v = self.hp.parse(vs)
        self.assertEqual(vs, str(v))
        self.assertEqual(v.posedit.pos.start.start.base, 90136803)
        self.assertEqual(v.posedit.pos.start.end.base, 90144453)
        self.assertEqual(v.posedit.pos.start.uncertain, True)
        self.assertEqual(v.posedit.pos.end.start.base, 90159675)
        self.assertEqual(v.posedit.pos.end.end.base, 90261231)
        self.assertEqual(v.posedit.pos.end.uncertain, True)
        self.assertEqual(type(v.posedit.edit).__name__, "Dup")

        vs2 = "NC_000009.11:g.(?_108337304)_(108337428_?)del"
        v2 = self.hp.parse(vs2)
        self.assertEqual(vs2, str(v2))
        self.assertEqual(v2.posedit.pos.start.start.base, None)
        self.assertEqual(v2.posedit.pos.start.uncertain, True)
        self.assertEqual(v2.posedit.pos.start.end.base, 108337304)
        self.assertEqual(v2.posedit.pos.end.start.base, 108337428)
        self.assertEqual(v2.posedit.pos.end.end.base, None)
        self.assertEqual(v2.posedit.pos.end.uncertain, True)
        self.assertEqual(type(v2.posedit.edit).__name__, "NARefAlt")

        # TODO add to parser test cases
        v3s = "NC_000005.9:g.(90136803_90159675)dup"
        v3 = self.hp.parse(v3s)

        # Interval itself is uncertain, but positions are certain
        self.assertEqual(v3s, str(v3))
        self.assertEqual(v3.posedit.pos.uncertain, True)
        self.assertEqual(v3.posedit.pos.start.uncertain, False)
        self.assertEqual(v3.posedit.pos.end.uncertain, False)
        self.assertEqual(v3.posedit.pos.start.base, 90136803)
        self.assertEqual(v3.posedit.pos.end.base, 90159675)


        v4s = "NC_000005.9:g.(90136803)_(90159675)dup"
        v4 = self.hp.parse(v4s)
        self.assertEqual(v4s, str(v4))
        self.assertEqual(v4.posedit.pos.uncertain, False)
        self.assertEqual(v4.posedit.pos.start.uncertain, True)
        self.assertEqual(str(v4.posedit.pos.start), "(90136803)")
        self.assertEqual(str(v4.posedit.pos.end), "(90159675)")
        self.assertEqual(v4.posedit.pos.end.uncertain, True)
        self.assertEqual(v4.posedit.pos.start.start.base, 90136803)
        self.assertEqual(v4.posedit.pos.end.start.base, 90159675)


        # Test that the start and end positions are not the same object
        v4.posedit.pos.start.end.base = 90136804
        self.assertTrue(v4.posedit.pos.start.start is not v4.posedit.pos.start.end)
        self.assertEqual(v4.posedit.pos.start.end.base, 90136804)
        self.assertEqual(v4.posedit.pos.start.start.base, 90136803)

        self.assertEqual("NC_000005.9:g.(90136803_90136804)_(90159675)dup", str(v4))


    def test_partial_uncertain_projection(self):
        data = [
            ("NC_000009.11:g.108337304_(108337428_?)del", False, True,
             "NM_001079802.1:n.207_(321+10_?)del",
             "NM_001079802.1:c.-10_(105+10_?)del"),

            ("NC_000009.11:g.(?_108337304)_108337428del", True, False,
             "NM_001079802.1:n.(?_207)_321+10del",
             "NM_001079802.1:c.(?_-10)_105+10del"),
        ]

        for hgvs_g, start_uncertain, stop_uncertain, hgvs_n, hgvs_c in data:
            var_g = self.hp.parse(hgvs_g)
            self.assertEqual(var_g.posedit.pos.start.uncertain, start_uncertain)
            self.assertEqual(var_g.posedit.pos.end.uncertain, stop_uncertain)
            self.assertEqual(hgvs_g, str(var_g))
            acc = hgvs_c.split(":")[0]
            var_n = self.vm.g_to_n(var_g, acc)
            self.assertEqual(hgvs_n, str(var_n))
            var_c = self.vm.g_to_c(var_g, acc)
            self.assertEqual(hgvs_c, str(var_c))

    def test_uncertain_projection_g_to_c_confidence(self):

        data = [
            ("NC_000005.9:g.(90136803_90144453)_(90159675_90261231)dup",
             "NM_032119.3:n.(?_17116-1)_(17952+1_?)dup",
             "NM_032119.3:c.(?_17020-1)_(17856+1_?)dup"),

            ("NC_000019.9:g.(11211022_11213339)_(11217364_11218067)dup",
             "NM_000527.5:n.(?_277-1)_(903+1_?)dup",
             "NM_000527.5:c.(?_191-1)_(817+1_?)dup"),

            ("NC_000009.11:g.(?_108337304)_(108337428_?)del",
             "NM_001079802.1:n.(?_207)_(321+10_?)del",
             "NM_001079802.1:c.(?_-10)_(105+10_?)del"),
        ]

        for hgvs_g, hgvs_n, hgvs_c in data:
            var_g = self.hp.parse(hgvs_g)
            self.assertEqual(hgvs_g, str(var_g))
            acc = hgvs_c.split(":")[0]
            var_n = self.vm.g_to_n(var_g, acc)
            self.assertEqual(hgvs_n, str(var_n))
            var_c = self.vm.g_to_c(var_g, acc)
            self.assertEqual(hgvs_c, str(var_c))


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
