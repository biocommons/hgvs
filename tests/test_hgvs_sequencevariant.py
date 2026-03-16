# -*- coding: utf-8 -*-
import pytest


import hgvs
import hgvs.parser
import hgvs.sequencevariant
import hgvs.normalizer


def test_gene_formatting(parser):
    v = parser.parse("NM_01234.5(BOGUS):c.65A>C")
    assert str(v) == "NM_01234.5(BOGUS):c.65A>C"
    v.gene = None
    assert str(v) == "NM_01234.5:c.65A>C"


@pytest.mark.vcr
@pytest.mark.quick
@pytest.mark.models
class Test_SequenceVariant:
    def test_SequenceVariant(self):
        var = hgvs.sequencevariant.SequenceVariant(
            ac="AC", type="B", posedit="1234DE>FG"
        )
        assert str(var) == "AC:B.1234DE>FG"

    def test_fill_ref(self, parser, hdp):
        hp = parser
        # fill reference for sequence variants
        var = hp.parse_hgvs_variant("NM_001166478.1:c.31_32del").fill_ref(hdp)
        assert var.format({"max_ref_length": None}) == "NM_001166478.1:c.31_32delTT"

        var = hp.parse_hgvs_variant("NM_001166478.1:c.31_32del2").fill_ref(hdp)
        assert var.format({"max_ref_length": None}) == "NM_001166478.1:c.31_32delTT"

        var = hp.parse_hgvs_variant("NM_001166478.1:c.2_7delinsTTTAGA").fill_ref(hdp)
        assert (
            var.format({"max_ref_length": None})
            == "NM_001166478.1:c.2_7delTGAAGAinsTTTAGA"
        )

        var = hp.parse_hgvs_variant("NM_001166478.1:c.35_36dup").fill_ref(hdp)
        assert var.format({"max_ref_length": None}) == "NM_001166478.1:c.35_36dupTC"

        var = hp.parse_hgvs_variant("NM_001166478.1:c.18_19insACT").fill_ref(hdp)
        assert var.format({"max_ref_length": None}) == "NM_001166478.1:c.18_19insACT"

        var = hp.parse_hgvs_variant("NM_001166478.1:c.31=").fill_ref(hdp)
        assert var.format({"max_ref_length": None}) == "NM_001166478.1:c.31T="

    def test_format(self, parser):
        hp = parser
        # Global default settings
        var = hp.parse_hgvs_variant("NP_001628.1:p.Gly528Arg")
        assert str(var) == "NP_001628.1:p.Gly528Arg"
        assert var.format() == "NP_001628.1:p.Gly528Arg"

        # Change global settings
        hgvs.global_config.formatting.p_3_letter = False
        assert str(var) == "NP_001628.1:p.G528R"

        # Custom settings
        hgvs.global_config.formatting.p_3_letter = True
        conf = {"p_3_letter": False}
        assert var.format(conf) == "NP_001628.1:p.G528R"

        var = hp.parse_hgvs_variant("NP_001628.1:p.Gly528Ter")
        conf = {"p_term_asterisk": True}
        assert var.format(conf) == "NP_001628.1:p.Gly528*"
        assert var.format() == "NP_001628.1:p.Gly528Ter"
        conf = {"p_3_letter": False}
        assert var.format(conf) == "NP_001628.1:p.G528*"
        assert var.format() == "NP_001628.1:p.Gly528Ter"

        # Remove reference sequence
        var = hp.parse_hgvs_variant("NM_001166478.1:c.31_32delTT")
        assert str(var) == "NM_001166478.1:c.31_32del"
        assert var.format(conf={"max_ref_length": 1}) == "NM_001166478.1:c.31_32del"
        assert var.format(conf={"max_ref_length": 2}) == "NM_001166478.1:c.31_32delTT"
        assert (
            var.format(conf={"max_ref_length": None}) == "NM_001166478.1:c.31_32delTT"
        )

        var = hp.parse_hgvs_variant("NM_001166478.1:c.31_32del2")
        assert str(var) == "NM_001166478.1:c.31_32del"
        assert var.format(conf={"max_ref_length": None}) == "NM_001166478.1:c.31_32del2"

        var = hp.parse_hgvs_variant("NM_001166478.1:c.31_32delTTinsAA")
        assert str(var) == "NM_001166478.1:c.31_32delinsAA"
        var = hp.parse_hgvs_variant("NM_001166478.1:c.35_36dupTC")
        assert str(var) == "NM_001166478.1:c.35_36dup"
        var = hp.parse_hgvs_variant("NM_001166478.1:c.31T=")
        assert str(var) == "NM_001166478.1:c.31="
        assert var.format(conf={"max_ref_length": None}) == "NM_001166478.1:c.31T="

    def test_uncertain(self, parser):
        hp = parser
        vs = "NC_000005.9:g.(90136803_90144453)_(90159675_90261231)dup"
        v = hp.parse(vs)
        assert vs == str(v)
        assert v.posedit.pos.start.start.base == 90136803
        assert v.posedit.pos.start.end.base == 90144453
        assert v.posedit.pos.start.uncertain is True
        assert v.posedit.pos.end.start.base == 90159675
        assert v.posedit.pos.end.end.base == 90261231
        assert v.posedit.pos.end.uncertain is True
        assert type(v.posedit.edit).__name__ == "Dup"

        vs2 = "NC_000009.11:g.(?_108337304)_(108337428_?)del"
        v2 = hp.parse(vs2)
        assert vs2 == str(v2)
        assert v2.posedit.pos.start.start.base is None
        assert v2.posedit.pos.start.uncertain is True
        assert v2.posedit.pos.start.end.base == 108337304
        assert v2.posedit.pos.end.start.base == 108337428
        assert v2.posedit.pos.end.end.base is None
        assert v2.posedit.pos.end.uncertain is True
        assert type(v2.posedit.edit).__name__ == "NARefAlt"

        # TODO add to parser test cases
        v3s = "NC_000005.9:g.(90136803_90159675)dup"
        v3 = hp.parse(v3s)

        # Interval itself is uncertain, but positions are certain
        assert v3s == str(v3)
        assert v3.posedit.pos.uncertain is True
        assert v3.posedit.pos.start.uncertain is False
        assert v3.posedit.pos.end.uncertain is False
        assert v3.posedit.pos.start.base == 90136803
        assert v3.posedit.pos.end.base == 90159675

        v4s = "NC_000005.9:g.(90136803)_(90159675)dup"
        v4 = hp.parse(v4s)
        assert v4s == str(v4)
        assert v4.posedit.pos.uncertain is False
        assert v4.posedit.pos.start.uncertain is True
        assert str(v4.posedit.pos.start) == "(90136803)"
        assert str(v4.posedit.pos.end) == "(90159675)"
        assert v4.posedit.pos.end.uncertain is True
        assert v4.posedit.pos.start.start.base == 90136803
        assert v4.posedit.pos.end.start.base == 90159675

        # Test that the start and end positions are not the same object
        v4.posedit.pos.start.end.base = 90136804
        assert v4.posedit.pos.start.start is not v4.posedit.pos.start.end
        assert v4.posedit.pos.start.end.base == 90136804
        assert v4.posedit.pos.start.start.base == 90136803

        assert str(v4) == "NC_000005.9:g.(90136803_90136804)_(90159675)dup"

    @pytest.mark.parametrize(
        "hgvs_g, start_uncertain, stop_uncertain, hgvs_n, hgvs_c",
        [
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
        ],
    )
    def test_partial_uncertain_projection(
        self, parser, vm, hdp, hgvs_g, start_uncertain, stop_uncertain, hgvs_n, hgvs_c
    ):
        hp = parser
        normalizer = hgvs.normalizer.Normalizer(hdp)
        """Test partial uncertain projection from genomic to cDNA to CDS coordinates."""
        var_g = hp.parse(hgvs_g)
        assert var_g.posedit.pos.start.uncertain == start_uncertain
        assert var_g.posedit.pos.end.uncertain == stop_uncertain
        assert hgvs_g == str(var_g)
        acc = hgvs_c.split(":")[0]
        var_n = vm.g_to_n(var_g, acc)
        assert hgvs_n == str(var_n)
        var_c = vm.g_to_c(var_g, acc)
        assert hgvs_c == str(var_c)

        var_g_norm = normalizer.normalize(var_g)
        assert hgvs_g == str(var_g_norm)

    def test_normalize_uncertain_hgvs_g(self, parser, vm, hdp):
        hp = parser
        normalizer = hgvs.normalizer.Normalizer(hdp)
        var_dup = "NC_000005.9:g.(90136802_90144453)_(90159674_90261231)dup"
        var_g = hp.parse(var_dup)
        var_g_norm = normalizer.normalize(var_g)
        assert str(var_g_norm) == var_dup
        assert var_g_norm.posedit.pos.start.uncertain
        assert var_g_norm.posedit.pos.end.uncertain
        assert var_g_norm.posedit.edit.type == "dup"

    @pytest.mark.parametrize(
        "hgvs_g, hgvs_n, hgvs_c",
        [
            (
                "NC_000005.9:g.(90136803_90144453)_(90159675_90261231)dup",
                "NM_032119.3:n.(17115+1_17116-1)_(17952+1_17953-1)dup",
                "NM_032119.3:c.(17019+1_17020-1)_(17856+1_17857-1)dup",
            ),
            (
                "NC_000019.9:g.(11211022_11213339)_(11217364_11218067)dup",
                "NM_000527.5:n.(276+1_277-1)_(903+1_904-1)dup",
                "NM_000527.5:c.(190+1_191-1)_(817+1_818-1)dup",
            ),
            (
                "NC_000009.11:g.(?_108337304)_(108337428_?)del",
                "NM_001079802.1:n.(?_207)_(321+10_?)del",
                "NM_001079802.1:c.(?_-10)_(105+10_?)del",
            ),
        ],
    )
    def test_uncertain_projection_g_to_c_confidence(
        self, parser, vm, hgvs_g, hgvs_n, hgvs_c
    ):
        hp = parser
        """Test uncertain projection from genomic to cDNA to CDS coordinates."""
        var_g = hp.parse(hgvs_g)
        assert hgvs_g == str(var_g)
        acc = hgvs_c.split(":")[0]
        var_n = vm.g_to_n(var_g, acc)
        assert hgvs_n == str(var_n)
        var_c = vm.g_to_c(var_g, acc)
        assert hgvs_c == str(var_c)

    @pytest.mark.parametrize(
        "clinvar_id, event_type, clinvar_hgvs_g, clinvar_hgvs_c",
        [
            (
                11692,
                "DEL",
                "NC_000023.11:g.(133661675_133661730)_(133661850_133661926)del",
                "NM_004484.3:c.(1293-76_1293)_(1413_1413+55)del",  # orig is this but it is wrong:"NM_004484.3:c.(1293_1293)-76_(1413_1413)del",
            ),
            (
                31562,
                "DEL",
                "NC_000007.14:g.(16361566_16366648)_(16369693_16391969)del",
                "NM_001101426.3:c.(534+14092_684+6399)_(684+9444_684+14526)del",  # orig is different: NM_001101426.3:c.(535_684)+6399_(535_684)+14526del
            ),
            (
                425669,
                "DEL",
                "NC_000002.12:g.(?_202376935)_(202377551_202464808)del",
                "NM_001204.6:c.(?_-540)_(76+1_77-1)del",
            ),
            (
                425698,
                "DEL",
                "NC_000002.12:g.(202377551_202464808)_(202559947_?)del",
                "NM_001204.6:c.(76+1_77-1)_(*1_?)del",
            ),
            (
                220591,
                "DEL",
                "NC_000017.11:g.(?_58709859)_(58734342_?)del",
                "NM_058216.2:c.(?_706)_(*120_?)del",  # orig seems different: NM_058216.2:c.706-?_*120del
            ),
            (
                251062,
                "DUP",
                "NC_000019.9:g.(11211022_11213339)_(11217364_11218067)dup",
                "NM_000527.5:c.(190+1_191-1)_(817+1_818-1)dup",
            ),
            (
                565301,
                "DUP",
                "NC_000001.11:g.(216073301_216078088)_(216327655_216364952)dup",
                "NM_206933.2:c.(784+1_785-1)_(5572+1_5573-1)dup",
            ),
            (
                254064,
                "DUP",
                "NC_000012.11:g.(?_133248801)_(133257865_?)dup",
                "NM_006231.3:c.(?_63)_(1794_?)dup",  # orig is different: NM_006231.3:c.63-?_1794+?dup1732
            ),
            (
                237630,
                "DUP",
                "NC_000017.10:g.(?_15133094)_(15164078_?)dup",
                "NM_000304.3:c.(?_-34)_(*1140_?)dup",  # orig is different NM_000304.3:c.-34-?_*1140dup1657
            ),
        ],
    )
    def test_clinvar_uncertain_ranges(
        self,
        parser,
        am38,
        vm,
        clinvar_id,
        event_type,
        clinvar_hgvs_g,
        clinvar_hgvs_c,
    ):
        """This is a unit test for the clinvar uncertain ranges described in
        issue #225: https://github.com/biocommons/hgvs/issues/225.

        This test goes in a loop -> uncertain hgvs_g -> hgvs_c -> hgvs_g.
        As part of this loop we lose information about the outer confidence interval, but we retain the inner confidence interval.
        That's why the hgvs_g_hgvs_g value contains only the inner interval.
        """
        hp = parser
        var_g = hp.parse(clinvar_hgvs_g)
        assert clinvar_hgvs_g == str(var_g)

        assert not var_g.posedit.pos.uncertain
        assert var_g.posedit.pos.start.uncertain

        chrom_ac = var_g.ac
        tx_ac = clinvar_hgvs_c.split(":")[0]

        var_c = vm.g_to_c(var_g, tx_ac)
        assert clinvar_hgvs_c == str(var_c)

        # the start/stop positions are uncertain, but the whole range is not:
        assert not var_c.posedit.pos.uncertain
        assert var_c.posedit.pos.start.uncertain
        assert var_c.posedit.pos.end.uncertain

        # var_n = vm.g_to_n(var_g, tx_ac)

        var_g_reverse = vm.c_to_g(var_c, chrom_ac)

        assert clinvar_hgvs_g == str(var_g_reverse)

        var_g_reverse_precise = vm.c_to_g(var_c, chrom_ac)
        assert clinvar_hgvs_g == str(var_g_reverse_precise)


if __name__ == "__main__":
    pytest.main()

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
