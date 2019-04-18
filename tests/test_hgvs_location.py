# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import unittest

import pytest

from hgvs.exceptions import HGVSError, HGVSUnsupportedOperationError
from hgvs.enums import Datum
import hgvs.location
import hgvs.parser


@pytest.mark.quick
@pytest.mark.models
class Test_SimplePosition(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hp = hgvs.parser.Parser()

    def test_success(self):
        self.assertEqual(str(hgvs.location.SimplePosition(5)), "5")
        self.assertEqual(str(hgvs.location.SimplePosition(5, uncertain=True)), "(5)")
        self.assertEqual(str(hgvs.location.SimplePosition(None)), "?")

    def test_failure(self):
        with self.assertRaises(AssertionError):
            self.assertEqual(hgvs.location.SimplePosition(-1), "SHOULD FAIL")

    def test_simple_subtraction(self):
        self.assertEqual(hgvs.location.SimplePosition(5) - hgvs.location.SimplePosition(3), 2)

    def test_simple_comparision(self):
        var = self.hp.parse_hgvs_variant("NC_000007.13:g.36561662_36561683del")
        self.assertFalse(var.posedit.pos.start == var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start < var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start <= var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start > var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start >= var.posedit.pos.end)

        var = self.hp.parse_hgvs_variant("NC_000007.13:g.36561662C>T")
        self.assertTrue(var.posedit.pos.start == var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start < var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start <= var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start > var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start >= var.posedit.pos.end)


@pytest.mark.quick
class Test_BaseOffsetPosition(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hp = hgvs.parser.Parser()

    def test_success(self):
        # r.5
        cdsp = hgvs.location.BaseOffsetPosition(5)
        self.assertEqual(cdsp.datum, Datum.SEQ_START)
        self.assertEqual(cdsp.base, 5)
        self.assertEqual(cdsp.offset, 0)
        self.assertEqual(str(cdsp), "5")
        self.assertFalse(cdsp.is_intronic)

        #r.5+6
        cdsp.offset = 6
        self.assertEqual(str(cdsp), "5+6")
        self.assertTrue(cdsp.is_intronic)

        #r.5+?
        cdsp.offset = None
        self.assertEqual(str(cdsp), "5+?")
        self.assertTrue(cdsp.is_intronic)

        #r.(5+?)
        cdsp.uncertain = True
        self.assertEqual(str(cdsp), "(5+?)")

        # c.*5
        cdsp = hgvs.location.BaseOffsetPosition(5, datum=Datum.CDS_END)
        self.assertEqual(cdsp.datum, Datum.CDS_END)
        self.assertEqual(cdsp.base, 5)
        self.assertEqual(cdsp.offset, 0)
        self.assertEqual(str(cdsp), "*5")

        cdsp.uncertain = True
        self.assertEqual(str(cdsp), "(*5)")

        cdsp.offset = 7
        self.assertEqual(str(cdsp), "(*5+7)")

    def test_baseoffset_subtraction(self):
        v30 = hgvs.location.BaseOffsetPosition(3, 0)
        v50 = hgvs.location.BaseOffsetPosition(5, 0)
        v52 = hgvs.location.BaseOffsetPosition(5, 2)
        v54 = hgvs.location.BaseOffsetPosition(5, 4)

        self.assertEqual(v50 - v30, 2)

        with self.assertRaises(HGVSError):
            _ = v54 - v30

    def test_baseoffset_comparision(self):
        var = self.hp.parse_hgvs_variant("NM_000030.2:c.669_680del")
        self.assertFalse(var.posedit.pos.start == var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start < var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start <= var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start > var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start >= var.posedit.pos.end)

        var = self.hp.parse_hgvs_variant("NM_000030.2:c.679_680+2del")
        self.assertFalse(var.posedit.pos.start == var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start < var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start <= var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start > var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start >= var.posedit.pos.end)

        var = self.hp.parse_hgvs_variant("NM_000030.2:c.-6_680+2del")
        self.assertFalse(var.posedit.pos.start == var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start < var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start <= var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start > var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start >= var.posedit.pos.end)

        var = self.hp.parse_hgvs_variant("NM_000030.2:c.680+2_680+10del")
        self.assertFalse(var.posedit.pos.start == var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start < var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start <= var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start > var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start >= var.posedit.pos.end)

        var = self.hp.parse_hgvs_variant("NM_000030.2:c.680+2_*82del")
        self.assertFalse(var.posedit.pos.start == var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start < var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start <= var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start > var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start >= var.posedit.pos.end)

        var = self.hp.parse_hgvs_variant("NM_000030.2:c.-12_*82del")
        self.assertFalse(var.posedit.pos.start == var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start < var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start <= var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start > var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start >= var.posedit.pos.end)

        var = self.hp.parse_hgvs_variant("NM_000030.2:c.680+2_681del")
        self.assertFalse(var.posedit.pos.start == var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start < var.posedit.pos.end)
        self.assertTrue(var.posedit.pos.start <= var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start > var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start >= var.posedit.pos.end)

        var = self.hp.parse_hgvs_variant("NM_000030.2:c.680+2_681-32del")
        with self.assertRaises(HGVSUnsupportedOperationError):
            var.posedit.pos.start < var.posedit.pos.end


@pytest.mark.quick
class Test_AAPosition(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hp = hgvs.parser.Parser()

    def test_AAPosition(self):
        ap = hgvs.location.AAPosition(15, "S")
        self.assertEqual(ap.pos, 15)
        self.assertEqual(str(ap), "Ser15")

    def test_aaposition_subtraction(self):
        l1 = hgvs.location.AAPosition(15, 'S')
        l2 = hgvs.location.AAPosition(20, 'S')
        self.assertEqual(l2 - l1, 5)

    def test_aaposition_comparision(self):
        var = self.hp.parse_hgvs_variant("NP_000042.3:p.His1082_Val1085delinsLeuHisGlnAla")
        self.assertTrue(var.posedit.pos.start < var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start > var.posedit.pos.end)

        var = self.hp.parse_hgvs_variant("NP_000042.3:p.His1082ArgfsTer2")
        self.assertFalse(var.posedit.pos.start < var.posedit.pos.end)
        self.assertFalse(var.posedit.pos.start > var.posedit.pos.end)


@pytest.mark.quick
class Test_Interval(unittest.TestCase):
    def test_Interval(self):
        ival = hgvs.location.Interval(
            hgvs.location.BaseOffsetPosition(base=12, offset=+34),
            hgvs.location.BaseOffsetPosition(base=56, offset=-78))
        self.assertEqual(ival.start.base, 12)
        self.assertEqual(ival.start.offset, 34)
        self.assertEqual(ival.end.base, 56)
        self.assertEqual(ival.end.offset, -78)
        self.assertEqual(str(ival), "12+34_56-78")

    def test_length(self):
        ival = hgvs.location.Interval(
            hgvs.location.BaseOffsetPosition(base=12, offset=0),
            hgvs.location.BaseOffsetPosition(base=50, offset=0))
        self.assertEqual(ival._length(), 39)


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
