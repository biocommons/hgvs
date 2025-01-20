# -*- coding: utf-8 -*-
import os
import unittest

import parameterized
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
        var = hgvs.sequencevariant.SequenceVariant(
            ac="AC", type="B", posedit="1234DE>FG"
        )
        self.assertEqual(str(var), "AC:B.1234DE>FG")

    def test_fill_ref(self):
        # fill reference for sequence variants
        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.31_32del").fill_ref(self.hdp)
        self.assertEqual(
            var.format({"max_ref_length": None}), "NM_001166478.1:c.31_32delTT"
        )

        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.31_32del2").fill_ref(
            self.hdp
        )
        self.assertEqual(
            var.format({"max_ref_length": None}), "NM_001166478.1:c.31_32delTT"
        )

        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.2_7delinsTTTAGA").fill_ref(
            self.hdp
        )
        self.assertEqual(
            var.format({"max_ref_length": None}),
            "NM_001166478.1:c.2_7delTGAAGAinsTTTAGA",
        )

        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.35_36dup").fill_ref(self.hdp)
        self.assertEqual(
            var.format({"max_ref_length": None}), "NM_001166478.1:c.35_36dupTC"
        )

        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.18_19insACT").fill_ref(
            self.hdp
        )
        self.assertEqual(
            var.format({"max_ref_length": None}), "NM_001166478.1:c.18_19insACT"
        )

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
        self.assertEqual(
            var.format(conf={"max_ref_length": 1}), "NM_001166478.1:c.31_32del"
        )
        self.assertEqual(
            var.format(conf={"max_ref_length": 2}), "NM_001166478.1:c.31_32delTT"
        )
        self.assertEqual(
            var.format(conf={"max_ref_length": None}), "NM_001166478.1:c.31_32delTT"
        )

        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.31_32del2")
        self.assertEqual(str(var), "NM_001166478.1:c.31_32del")
        self.assertEqual(
            var.format(conf={"max_ref_length": None}), "NM_001166478.1:c.31_32del2"
        )

        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.31_32delTTinsAA")
        self.assertEqual(str(var), "NM_001166478.1:c.31_32delinsAA")
        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.35_36dupTC")
        self.assertEqual(str(var), "NM_001166478.1:c.35_36dup")
        var = self.hp.parse_hgvs_variant("NM_001166478.1:c.31T=")
        self.assertEqual(str(var), "NM_001166478.1:c.31=")
        self.assertEqual(
            var.format(conf={"max_ref_length": None}), "NM_001166478.1:c.31T="
        )

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
            (
                "NC_000009.11:g.108337304_(108337428_?)del",
                False,
                True,
                "NM_001079802.1:n.207_(321+10_?)del",
                "NM_001079802.1:c.-10_(105+10_?)del",
            ),
            (
                "NC_000009.11:g.(?_108337304)_108337428del",
                True,
                False,
                "NM_001079802.1:n.(?_207)_321+10del",
                "NM_001079802.1:c.(?_-10)_105+10del",
            ),
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
            (
                "NC_000005.9:g.(90136803_90144453)_(90159675_90261231)dup",
                "NM_032119.3:n.(?_17116-1)_(17952+1_?)dup",
                "NM_032119.3:c.(?_17020-1)_(17856+1_?)dup",
            ),
            (
                "NC_000019.9:g.(11211022_11213339)_(11217364_11218067)dup",
                "NM_000527.5:n.(?_277-1)_(903+1_?)dup",
                "NM_000527.5:c.(?_191-1)_(817+1_?)dup",
            ),
            (
                "NC_000009.11:g.(?_108337304)_(108337428_?)del",
                "NM_001079802.1:n.(?_207)_(321+10_?)del",
                "NM_001079802.1:c.(?_-10)_(105+10_?)del",
            ),
        ]

        for hgvs_g, hgvs_n, hgvs_c in data:
            var_g = self.hp.parse(hgvs_g)
            self.assertEqual(hgvs_g, str(var_g))
            acc = hgvs_c.split(":")[0]
            var_n = self.vm.g_to_n(var_g, acc)
            self.assertEqual(hgvs_n, str(var_n))
            var_c = self.vm.g_to_c(var_g, acc)
            self.assertEqual(hgvs_c, str(var_c))

    @parameterized.parameterized.expand(
        [
            (
                11692,
                "DEL",
                "NC_000023.11:g.(133661675_133661730)_(133661850_133661926)del",
                "NC_000023.11:g.(?_133661730)_(133661850_?)del",
                "NM_004484.3:c.(1293_1293)-76_(1413_1413)del",
                "NM_004484.3:c.(?_1293)_(1413_?)del",
            ),
            (
                31562,
                "DEL",
                "NC_000007.14:g.(16361566_16366648)_(16369693_16391969)del",
                "NC_000007.14:g.(?_16366648)_(16369693_?)del",
                "NM_001101426.3:c.(535_684)+6399_(535_684)+14526del",
                "NM_001101426.3:c.(?_684+6399)_(684+9444_?)del",
            ),
            (
                425669,
                "DEL",
                "NC_000002.12:g.(?_202376935)_(202377551_202464808)del",
                "NC_000002.12:g.(?_202376935)_(202377551_?)del",
                "NM_001204.6:c.(?_-540)_(76+1_77-1)del",
                "NM_001204.6:c.(?_-540)_(76+1_?)del",
            ),
            (
                425698,
                "DEL",
                "NC_000002.12:g.(202377551_202464808)_(202559947_?)del",
                "NC_000002.12:g.(?_202464808)_(202559947_?)del",
                "NM_001204.6:c.(76+1_77-1)_(*1_?)del",
                "NM_001204.6:c.(?_77-1)_(*1_?)del",
            ),
            (
                220591,
                "DEL",
                "NC_000017.11:g.(?_58709859)_(58734342_?)del",
                "NC_000017.11:g.(?_58709859)_(58734342_?)del",
                "NM_058216.2:c.706-?_*120del",
                "NM_058216.2:c.(?_706)_(*120_?)del",
            ),
            (
                251062,
                "DUP",
                "NC_000019.9:g.(11211022_11213339)_(11217364_11218067)dup",
                "NC_000019.9:g.(?_11213339)_(11217364_?)dup",
                "NM_000527.5:c.(190+1_191-1)_(817+1_818-1)dup",
                "NM_000527.5:c.(?_191-1)_(817+1_?)dup",
            ),
            (
                565301,
                "DUP",
                "NC_000001.11:g.(216073301_216078088)_(216327655_216364952)dup",
                "NC_000001.11:g.(?_216078088)_(216327655_?)dup",
                "NM_206933.2:c.(784+1_785-1)_(5572+1_5573-1)dup",
                "NM_206933.2:c.(?_785-1)_(5572+1_?)dup",
            ),
            (
                254064,
                "DUP",
                "NC_000012.11:g.(?_133248801)_(133257865_?)dup",
                "NC_000012.11:g.(?_133248801)_(133257865_?)dup",
                "NM_006231.3:c.63-?_1794+?dup1732",
                "NM_006231.3:c.(?_63)_(1794_?)dup",
            ),
            (
                237630,
                "DUP",
                "NC_000017.10:g.(?_15133094)_(15164078_?)dup",
                "NC_000017.10:g.(?_15133094)_(15164078_?)dup",
                "NM_000304.3:c.-34-?_*1140dup1657",
                "NM_000304.3:c.(?_-34)_(*1140_?)dup",
            ),
        ]
    )
    def test_clinvar_uncertain_ranges(
        self,
        clinvar_id,
        event_type,
        clinvar_hgvs_g,
        hgvs_hgvs_g,
        clinvar_hgvs_c,
        hgvs_hgvs_c,
    ):
        """This is a unit test for the clinvar uncertain ranges described in
        issue #225: https://github.com/biocommons/hgvs/issues/225.

        This test goes in a loop -> uncertain hgvs_g -> hgvs_c -> hgvs_g.
        As part of this loop we lose information about the outer confidence interval, but we retain the inner confidence internval.
        That's why the hgvs_g_hgvs_g value contains only the inner interval.
        """

        var_g = self.hp.parse(clinvar_hgvs_g)
        self.assertEqual(clinvar_hgvs_g, str(var_g))

        assert not var_g.posedit.pos.uncertain
        assert var_g.posedit.pos.start.uncertain

        chrom_ac = var_g.ac
        tx_ac = clinvar_hgvs_c.split(":")[0]
        var_c = self.vm.g_to_c(var_g, tx_ac)
        self.assertEqual(hgvs_hgvs_c, str(var_c))

        # the start/stop positions are uncertain, but the whole range is not:
        assert not var_c.posedit.pos.uncertain
        assert var_c.posedit.pos.start.uncertain
        assert var_c.posedit.pos.end.uncertain

        var_g_reverse = self.vm.c_to_g(var_c, chrom_ac)
        self.assertEqual(hgvs_hgvs_g, str(var_g_reverse))


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
