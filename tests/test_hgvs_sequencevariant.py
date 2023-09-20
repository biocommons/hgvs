# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

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
    def test_SequenceVariant(self):
        var = hgvs.sequencevariant.SequenceVariant(ac="AC", type="B", posedit="1234DE>FG")
        self.assertEqual(str(var), "AC:B.1234DE>FG")

    def test_fill_ref(self):
        hp = hgvs.parser.Parser()
        hdp = hgvs.dataproviders.uta.connect(
            mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE
        )

        # fill reference for sequence variants
        var = hp.parse_hgvs_variant("NM_001166478.1:c.31_32del").fill_ref(hdp)
        self.assertEqual(var.format({"max_ref_length": None}), "NM_001166478.1:c.31_32delTT")

        var = hp.parse_hgvs_variant("NM_001166478.1:c.31_32del2").fill_ref(hdp)
        self.assertEqual(var.format({"max_ref_length": None}), "NM_001166478.1:c.31_32delTT")

        var = hp.parse_hgvs_variant("NM_001166478.1:c.2_7delinsTTTAGA").fill_ref(hdp)
        self.assertEqual(
            var.format({"max_ref_length": None}), "NM_001166478.1:c.2_7delTGAAGAinsTTTAGA"
        )

        var = hp.parse_hgvs_variant("NM_001166478.1:c.35_36dup").fill_ref(hdp)
        self.assertEqual(var.format({"max_ref_length": None}), "NM_001166478.1:c.35_36dupTC")

        var = hp.parse_hgvs_variant("NM_001166478.1:c.18_19insACT").fill_ref(hdp)
        self.assertEqual(var.format({"max_ref_length": None}), "NM_001166478.1:c.18_19insACT")

        var = hp.parse_hgvs_variant("NM_001166478.1:c.31=").fill_ref(hdp)
        self.assertEqual(var.format({"max_ref_length": None}), "NM_001166478.1:c.31T=")

    def test_format(self):
        hp = hgvs.parser.Parser()

        # Global default settings
        var = hp.parse_hgvs_variant("NP_001628.1:p.Gly528Arg")
        self.assertEqual(str(var), "NP_001628.1:p.Gly528Arg")
        self.assertEqual(var.format(), "NP_001628.1:p.Gly528Arg")

        # Change global settings
        hgvs.global_config.formatting.p_3_letter = False
        self.assertEqual(str(var), "NP_001628.1:p.G528R")

        # Custom settings
        hgvs.global_config.formatting.p_3_letter = True
        conf = {"p_3_letter": False}
        self.assertEqual(var.format(conf), "NP_001628.1:p.G528R")

        var = hp.parse_hgvs_variant("NP_001628.1:p.Gly528Ter")
        conf = {"p_term_asterisk": True}
        self.assertEqual(var.format(conf), "NP_001628.1:p.Gly528*")
        self.assertEqual(var.format(), "NP_001628.1:p.Gly528Ter")
        conf = {"p_3_letter": False}
        self.assertEqual(var.format(conf), "NP_001628.1:p.G528*")
        self.assertEqual(var.format(), "NP_001628.1:p.Gly528Ter")

        # Remove reference sequence
        var = hp.parse_hgvs_variant("NM_001166478.1:c.31_32delTT")
        self.assertEqual(str(var), "NM_001166478.1:c.31_32del")
        self.assertEqual(var.format(conf={"max_ref_length": 1}), "NM_001166478.1:c.31_32del")
        self.assertEqual(var.format(conf={"max_ref_length": 2}), "NM_001166478.1:c.31_32delTT")
        self.assertEqual(var.format(conf={"max_ref_length": None}), "NM_001166478.1:c.31_32delTT")

        var = hp.parse_hgvs_variant("NM_001166478.1:c.31_32del2")
        self.assertEqual(str(var), "NM_001166478.1:c.31_32del")
        self.assertEqual(var.format(conf={"max_ref_length": None}), "NM_001166478.1:c.31_32del2")

        var = hp.parse_hgvs_variant("NM_001166478.1:c.31_32delTTinsAA")
        self.assertEqual(str(var), "NM_001166478.1:c.31_32delinsAA")
        var = hp.parse_hgvs_variant("NM_001166478.1:c.35_36dupTC")
        self.assertEqual(str(var), "NM_001166478.1:c.35_36dup")
        var = hp.parse_hgvs_variant("NM_001166478.1:c.31T=")
        self.assertEqual(str(var), "NM_001166478.1:c.31=")
        self.assertEqual(var.format(conf={"max_ref_length": None}), "NM_001166478.1:c.31T=")

    def test_uncertain(self):
        hp = hgvs.parser.Parser()

        vs = "NC_000005.9:g.(90136803_90144453)_(90159675_90261231)dup"
        v = hp.parse(vs)
        self.assertEqual(vs, str(v))
        self.assertEqual(v.posedit.pos.start.start.base, 90136803)
        self.assertEqual(v.posedit.pos.start.end.base, 90144453)
        self.assertEqual(v.posedit.pos.end.start.base, 90159675)
        self.assertEqual(v.posedit.pos.end.end.base, 90261231)
        self.assertEqual(type(v.posedit.edit).__name__, "Dup")

        vs2 = "NC_000009.11:g.(?_108337304)_(108337428_?)del"
        v2 = hp.parse(vs2)
        self.assertEqual(vs2, str(v2))
        self.assertEqual(v2.posedit.pos.start.start.base, None)
        self.assertEqual(v2.posedit.pos.start.uncertain, True)
        self.assertEqual(v2.posedit.pos.start.end.base, 108337304)
        self.assertEqual(v2.posedit.pos.end.start.base, 108337428)
        self.assertEqual(v2.posedit.pos.end.end.base, None)
        self.assertEqual(v2.posedit.pos.end.uncertain, True)
        self.assertEqual(type(v2.posedit.edit).__name__, "NARefAlt")



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
