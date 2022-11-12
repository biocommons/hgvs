# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import unittest

import pytest
from support import CACHE

import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper
from hgvs.exceptions import HGVSInvalidVariantError


@pytest.mark.quick
class Test_VariantMapper(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hdp = hgvs.dataproviders.uta.connect(mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE)
        cls.am = hgvs.assemblymapper.AssemblyMapper(cls.hdp)
        cls.am37 = hgvs.assemblymapper.AssemblyMapper(cls.hdp, assembly_name="GRCh37")
        cls.hp = hgvs.parser.Parser()

    def test_VariantMapper_quick(self):
        # From garcia.tsv:
        hgvs_g = "NC_000007.13:g.36561662C>T"
        hgvs_c = "NM_001637.3:c.1582G>A"
        hgvs_p = "NP_001628.1:p.(Gly528Arg)"

        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_c = self.am.g_to_c(var_g, "NM_001637.3")
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_c), hgvs_c)
        self.assertEqual(str(var_p), hgvs_p)

        with self.assertRaises(HGVSInvalidVariantError):
            self.am.c_to_p(self.hp.parse_hgvs_variant("NM_000059.3:c.7790delAAG"))

        hgvs_c = "NM_000059.3:c.7791A>G"
        hgvs_p = "NP_000050.2:p.(Lys2597=)"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_p), hgvs_p)

        hgvs_c = "NM_000302.3:c.1594_1596del"
        hgvs_p = "NP_000293.2:p.(Glu532del)"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_p), hgvs_p)

        hgvs_c = "NM_000090.3:c.2490_2516del"
        hgvs_p = "NP_000081.1:p.(Glu832_Gly840del)"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_p), hgvs_p)

        hgvs_c = "NM_001292004.1:c.376="
        hgvs_g = "NC_000005.10:g.74715659="

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_g = self.am.c_to_g(var_c)

        self.assertEqual(str(var_g), hgvs_g)

        hgvs_c = "NM_000116.4:c.-8_-3inv6"
        hgvs_g = "NC_000023.11:g.154411836_154411841inv"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_g = self.am.c_to_g(var_c)

        self.assertEqual(str(var_g), hgvs_g)

        hgvs_c = "NM_000348.3:c.89_91inv3"
        hgvs_g = "NC_000002.12:g.31580810_31580812inv"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_g = self.am.c_to_g(var_c)

        self.assertEqual(str(var_g), hgvs_g)

        hgvs_c = "NM_001637.3:c.1582_1583inv"
        hgvs_p = "NP_001628.1:p.(Gly528Pro)"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_p), hgvs_p)

        hgvs_c = "NM_025137.3:c.-20_*20inv"
        hgvs_p = "NP_079413.3:p.?"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_p), hgvs_p)

    def test_projection_at_alignment_discrepancy(self):
        hgvs_g = "NC_000019.10:g.50378563_50378564insTG"
        hgvs_n = "NM_007121.5:n.796_798delinsTG"

        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_n = self.am.g_to_n(var_g, "NM_007121.5")

        self.assertEqual(str(var_n), hgvs_n)

        hgvs_g = "NC_000007.14:g.149779575delC"
        hgvs_n = "NM_198455.2:n.1115_1116insAG"

        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_n = self.am.g_to_n(var_g, "NM_198455.2")

        self.assertEqual(str(var_n), hgvs_n)

        # issue-353
        hgvs_g = "NC_000012.11:g.122064775C>T"
        hgvs_c = "NM_032790.3:c.127_128insTGCCAC"

        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_c = self.am.g_to_c(var_g, "NM_032790.3")

        self.assertEqual(str(var_c), hgvs_c)

        # issue-461
        hgvs_g = "NC_000002.11:g.73675227_73675228insCTC"
        hgvs_c = "NM_015120.4:c.1574_1576="

        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_c = self.am37.g_to_c(var_g, "NM_015120.4")

        self.assertEqual(str(var_c), hgvs_c)

        # issue-259
        hgvs_c = "NM_000116.4:c.-120_-119insT"
        hgvs_g = "NC_000023.10:g.153640061="

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_g = self.am37.c_to_g(var_c)

        self.assertEqual(str(var_g), hgvs_g)

        hgvs_c = "NM_000348.3:c.88del"
        hgvs_g = "NC_000002.11:g.31805882_31805883="

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_g = self.am37.c_to_g(var_c)

        self.assertEqual(str(var_g), hgvs_g)

    def test_c_to_p_with_stop_gain(self):
        # issue-474
        hgvs_c = "NM_080877.2:c.1733_1735delinsTTT"
        hgvs_p = "NP_543153.1:p.(Pro578_Lys579delinsLeuTer)"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_p), hgvs_p)

        # issue-492
        hgvs_c = "NM_001034853.1:c.2847_2848delAGinsCT"
        hgvs_p = "NP_001030025.1:p.(Glu949_Glu950delinsAspTer)"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_p), hgvs_p)

        hgvs_c = "NM_001034853.1:c.2847_2848inv"
        hgvs_p = "NP_001030025.1:p.(Glu949_Glu950delinsAspTer)"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_p), hgvs_p)

        hgvs_c = "NM_080877.2:c.1735A>T"
        hgvs_p = "NP_543153.1:p.(Lys579Ter)"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_p), hgvs_p)

        hgvs_c = "NM_080877.2:c.1795_*3delinsTAG"
        hgvs_p = "NP_543153.1:p.(Leu599Ter)"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_p), hgvs_p)


class Test_RefReplacement(unittest.TestCase):
    test_cases = [
        # These casese attempt to test reference update in four dimensions:
        # - variant type: n, c, g
        # - major mapping paths: c<->n, c<->g, n<->g
        # - variant class: sub, del, ins, delins, dup
        # - strand: +/-
        # ADRB2    │ NM_000024.5 │  239 │ 1481 │ NC_000005.9  │  1 │ 148206155,148208197 | 284=1X32=1X1724=
        # cseq = hdp.fetch_seq("NM_000024.5")
        # gseq = hdp.fetch_seq("NC_000005.9",148206155,148208197)
        # cseq[280:290] = "CAATAGAAGC"
        # gseq[280:290] = "CAATGGAAGC"
        #                      ^ @ n.285
        # These variants are in and around the first sub:
        {"c": "NM_000024.5:c.42C>N", "g": "NC_000005.9:g.148206436C>N", "n": "NM_000024.5:n.281C>N"},
        {"c": "NM_000024.5:c.43A>N", "g": "NC_000005.9:g.148206437A>N", "n": "NM_000024.5:n.282A>N"},
        {"c": "NM_000024.5:c.44A>N", "g": "NC_000005.9:g.148206438A>N", "n": "NM_000024.5:n.283A>N"},
        {"c": "NM_000024.5:c.45T>N", "g": "NC_000005.9:g.148206439T>N", "n": "NM_000024.5:n.284T>N"},
        {"c": "NM_000024.5:c.46A>N", "g": "NC_000005.9:g.148206440G>N", "n": "NM_000024.5:n.285A>N"},  # ref repl
        {"c": "NM_000024.5:c.47G>N", "g": "NC_000005.9:g.148206441G>N", "n": "NM_000024.5:n.286G>N"},
        {"c": "NM_000024.5:c.48A>N", "g": "NC_000005.9:g.148206442A>N", "n": "NM_000024.5:n.287A>N"},
        {"c": "NM_000024.5:c.49A>N", "g": "NC_000005.9:g.148206443A>N", "n": "NM_000024.5:n.288A>N"},
        {"c": "NM_000024.5:c.50G>N", "g": "NC_000005.9:g.148206444G>N", "n": "NM_000024.5:n.289G>N"},
        {"c": "NM_000024.5:c.51C>N", "g": "NC_000005.9:g.148206445C>N", "n": "NM_000024.5:n.290C>N"},
        # ins, del, delins, dup:
        {
            "c": "NM_000024.5:c.46_47insNN",
            "g": "NC_000005.9:g.148206440_148206441insNN",
            "n": "NM_000024.5:n.285_286insNN",
        },
        {
            "c": "NM_000024.5:c.45_47delTAG",
            "g": "NC_000005.9:g.148206439_148206441delTGG",
            "n": "NM_000024.5:n.284_286delTAG",
        },
        {
            "c": "NM_000024.5:c.45_47delTAGinsNNNN",
            "g": "NC_000005.9:g.148206439_148206441delTGGinsNNNN",
            "n": "NM_000024.5:n.284_286delTAGinsNNNN",
        },
        {
            "c": "NM_000024.5:c.45_47delTAGinsNNNN",
            "g": "NC_000005.9:g.148206439_148206441delTGGinsNNNN",
            "n": "NM_000024.5:n.284_286delTAGinsNNNN",
        },
        {"c": "NM_000024.5:c.46dupA", "g": "NC_000005.9:g.148206440dupG", "n": "NM_000024.5:n.285dupA"},
        # IFNA16   │ NM_002173.2 │    6 │  576 │ NC_000009.11 │ -1 │  21216371, 21217310 | 691=2X246=
        # cseq = hdp.fetch_seq("NM_002173.2")
        # gseq = reverse_complement(hdp.fetch_seq("NC_000009.11",21216371,21217310))
        # cseq[685:695] = "AAATTTCAAA"
        # gseq[685:695] = "AAATTTTCAA"
        #                        ^^ @ n.692_693
        # These variants are in and around the 2X substitution
        {"c": "NM_002173.2:c.*110A>N", "g": "NC_000009.11:g.21216625T>N", "n": "NM_002173.2:n.686A>N"},
        {"c": "NM_002173.2:c.*111A>N", "g": "NC_000009.11:g.21216624T>N", "n": "NM_002173.2:n.687A>N"},
        {"c": "NM_002173.2:c.*112A>N", "g": "NC_000009.11:g.21216623T>N", "n": "NM_002173.2:n.688A>N"},
        {"c": "NM_002173.2:c.*113T>N", "g": "NC_000009.11:g.21216622A>N", "n": "NM_002173.2:n.689T>N"},
        {"c": "NM_002173.2:c.*114T>N", "g": "NC_000009.11:g.21216621A>N", "n": "NM_002173.2:n.690T>N"},
        {"c": "NM_002173.2:c.*115T>N", "g": "NC_000009.11:g.21216620A>N", "n": "NM_002173.2:n.691T>N"},
        {"c": "NM_002173.2:c.*116C>N", "g": "NC_000009.11:g.21216619A>N", "n": "NM_002173.2:n.692C>N"},  # ref repl
        {"c": "NM_002173.2:c.*117A>N", "g": "NC_000009.11:g.21216618G>N", "n": "NM_002173.2:n.693A>N"},  # ref repl
        {"c": "NM_002173.2:c.*118A>N", "g": "NC_000009.11:g.21216617T>N", "n": "NM_002173.2:n.694A>N"},
        {"c": "NM_002173.2:c.*119A>N", "g": "NC_000009.11:g.21216616T>N", "n": "NM_002173.2:n.695A>N"},
        # ins, del, delins, dup:
        {
            "c": "NM_002173.2:c.*115_*117insNN",
            "g": "NC_000009.11:g.21216618_21216620insNN",
            "n": "NM_002173.2:n.691_693insNN",
        },
        {
            "c": "NM_002173.2:c.*114_*117delTTCA",
            "g": "NC_000009.11:g.21216618_21216621delGAAA",
            "n": "NM_002173.2:n.690_693delTTCA",
        },
        {
            "c": "NM_002173.2:c.*115_*117delTCAinsNN",
            "g": "NC_000009.11:g.21216618_21216620delGAAinsNN",
            "n": "NM_002173.2:n.691_693delTCAinsNN",
        },
        {
            "c": "NM_002173.2:c.*115_*117delTCAinsNN",
            "g": "NC_000009.11:g.21216618_21216620delGAAinsNN",
            "n": "NM_002173.2:n.691_693delTCAinsNN",
        },
        {
            "c": "NM_002173.2:c.*115_*117dupTCA",
            "g": "NC_000009.11:g.21216618_21216620dupGAA",
            "n": "NM_002173.2:n.691_693dupTCA",
        },
        # genomic variation within an indel discrepancy
        # NM_032790.3  c        108 >
        # NM_032790.3  n        301 > CGGGGAGCCCCCGGGGGCC------CCGCCACCGCCGCCGT
        #                             |||||||||||||||||||------||||||||||||||||
        # NC_000012.11 g  122064755 > CGGGGAGCCCCCGGGGGCCCCGCCACCGCCACCGCCGCCGT
        #                                              **********
        {
            "c": "NM_032790.3:c.125_128delCCCC",
            "g": "NC_000012.11:g.122064772_122064781delCCCCGCCACC",
            "n": "NM_032790.3:n.318_321delCCCC",
        },
    ]

    @classmethod
    def setUpClass(cls):
        def _parse_rec(rec):
            rec["pv"] = {x: cls.hp.parse_hgvs_variant(rec[x]) for x in "cgn"}
            return rec

        cls.hdp = hgvs.dataproviders.uta.connect(mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE)
        cls.am = hgvs.assemblymapper.AssemblyMapper(
            cls.hdp, replace_reference=True, assembly_name="GRCh37", alt_aln_method="splign"
        )
        cls.hp = hgvs.parser.Parser()
        cls.tests = [_parse_rec(rec) for rec in cls.test_cases]

    def test_replace_reference_sequence(self):
        """AssemblyMapper: Replace invalid reference sequence"""

        for rec in self.tests:
            for x in "cgn":
                pv = rec["pv"][x]
                if pv.posedit.edit.ref:
                    # replace ref with junk
                    pv.posedit.edit.ref = "NNNNNN"
                self.am._replace_reference(pv)
                self.assertEqual(rec[x], pv.format(conf={"max_ref_length": None}))


@pytest.mark.quick
class Test_AssemblyMapper(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        hdp = hgvs.dataproviders.uta.connect(mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE)
        cls.hp = hgvs.parser.Parser()
        cls.am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37", alt_aln_method="splign")

    def _test_mapping(self, hgvs_set):
        """given list of variant strings, test all valid combinations of
        g<->n<->c<->p mappings

        """
        parsed_variants = [(hv, self.hp.parse_hgvs_variant(hv)) for hv in hgvs_set]
        hgvs = {v.type: hv for hv, v in parsed_variants}
        pvs = {v.type: v for hv, v in parsed_variants}

        if "g" in pvs and "c" in pvs:
            self.assertEqual(hgvs["g"], str(self.am.c_to_g(pvs["c"])))
            self.assertEqual(hgvs["c"], str(self.am.g_to_c(pvs["g"], pvs["c"].ac)))
        if "g" in pvs and "n" in pvs:
            self.assertEqual(hgvs["g"], str(self.am.n_to_g(pvs["n"])))
            self.assertEqual(hgvs["n"], str(self.am.g_to_n(pvs["g"], pvs["n"].ac)))
        if "c" in pvs and "n" in pvs:
            self.assertEqual(hgvs["n"], str(self.am.c_to_n(pvs["c"])))
            self.assertEqual(hgvs["c"], str(self.am.n_to_c(pvs["n"])))
        if "c" in pvs and "p" in pvs:
            self.assertEqual(hgvs["p"], str(self.am.c_to_p(pvs["c"])))
            self.assertEqual(hgvs["p"], str(self.am.t_to_p(pvs["c"])))

    def test_SNV(self):
        """AssemblyMapper: smoketest with SNVs"""
        hgvs_set = [
            "NC_000007.13:g.36561662C>T",
            "NM_001637.3:c.1582G>A",
            "NM_001637.3:n.1983G>A",
            "NP_001628.1:p.(Gly528Arg)",
        ]
        self._test_mapping(hgvs_set)

    def test_intronic(self):
        """AssemblyMapper: smoketest with intronic SNVs"""
        hgvs_set = [
            "NC_000010.10:g.89711873A>C",
            "NM_000314.4:c.493-2A>C",
            "NM_000314.4:n.1524-2A>C",
            "NP_000305.3:p.?",
        ]
        self._test_mapping(hgvs_set)

    def test_t_to_p(self):
        assert "non-coding" == str(self.am.t_to_p(self.hp.parse("NR_027676.1:n.3980del")))
        assert "NP_000050.2:p.(Lys2597=)" == str(self.am.t_to_p(self.hp.parse("NM_000059.3:c.7791A>G")))


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
