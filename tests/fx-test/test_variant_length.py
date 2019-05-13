# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import os

import unittest

import pytest

from hgvs.exceptions import HGVSUnsupportedOperationError
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper
from support import CACHE


class Test_VariantLengths(unittest.TestCase):
    """test length_change method for all variant types and cases"""

    success_tests = [
    # class AAExt(Edit),
        ("NP_000040.1:p.Ter314Trpext*45", 45),
        ("NP_000040.1:p.Met1ext-12", 12),

    # class AAFs(Edit):
    # UNSUPPORTED: ("NP_000001.1:p.Gly25Trpfs*10", ),   # requires seq fetch for length :-(

    # class AARefAlt(Edit):
        ("NP_006149.2:p.(Glu396delinsLysLys)", 1),    # NM_006158.3:c.1186delGinsAAAA
        ("NP_006149.2:p.(Glu397del)", -1),    # NM_006158.3:n.1288_1290delGAG
        ("NP_006149.2:p.(Glu397dup)", 1),    # NM_006158.3:n.1288_1290delGAGinsGAGGAG

    # class AASub(AARefAlt):
        ("NP_006149.2:p.(Glu396Lys)", 0),    # NM_006158.3:c.1186delGinsA

    # class Conv(Edit),
    # UNSUPPORTED

    # class Dup(Edit):
        ("NM_006158.3:n.1291_1293dup", 3),
        ("NM_006158.3:n.1291_1293dupGAG", 3),

    # class Inv(Edit):
        ("NM_006158.3:n.1291_1293inv", 0),

    # class NACopy(Edit):
        ("NC_000014.8:g.88401076_88459508copy4", (88459508 - 88401076 + 1) * 4),

    # class    NARefAlt(Edit):
        ("NM_000314.4:c.706G>T", 0),
        ("NM_000314.4:c.706delG", -1),
        ("NM_000314.4:c.706delGinsT", 0),
        ("NM_000314.4:c.706_708del", -3),
        ("NM_000314.4:c.706_706delG", -1),
        ("NM_000314.4:c.706_708delGAC", -3),
        ("NM_000314.4:c.706_708delGACinsT", -2),
        ("NM_000314.4:c.706_708delGACinsTT", -1),
        ("NM_000314.4:c.706_708delGACinsTTG", 0),
        ("NM_000314.4:c.706_708delGACinsTTGT", +1),
        ("NM_000314.4:c.706_708del3insTTGT", +1),
        ("NM_000314.4:c.706_708del3", -3),
        ("NM_004655.3:c.-12_8delAGCTCCCTCACCATGAGTAG", -20),
        ("NM_000314.4:c.706_707insT", +1),
        ("NM_000314.4:c.706_707insTT", +2),
    #See #319:: ("NM_000249.3:c.1897-3del", ),
    #See #319:: ("NM_000249.3:c.1897-3_1897-1del", ),
    #See #319:: ("NM_000249.3:c.1897-3_1897-1delinsACGT", ),

    # identity variants:
        ("NM_000314.4:c.1=", 0),
        ("NM_000314.4:c.1C=", 0),
        ("NM_000314.4:c.1_2CC=", 0),
        ("NM_000314.4:c.1_2delinsCC", 0),
    ]

    error_tests = [
        "NM_003777.3:c.13552_*36del57",
    ]

    @classmethod
    def setUpClass(cls):
        cls.hdp = hgvs.dataproviders.uta.connect(
            mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE)
        cls.vm = hgvs.variantmapper.VariantMapper(cls.hdp)
        cls.hp = hgvs.parser.Parser()

    def test_success_cases(self):
        "posedit length_change: test supported variant types"
        for tv, tl in self.success_tests:
            v = self.hp.parse_hgvs_variant(tv)
            self.assertEqual(v.posedit.length_change(), tl)

    def test_error_cases_w_exceptions(self):
        "posedit length_change: raise exception on error"
        for tv in self.error_tests:
            v = self.hp.parse_hgvs_variant(tv)
            with self.assertRaises(HGVSUnsupportedOperationError):
                _ = v.posedit.length_change()

    def test_error_cases_w_error_values(self):
        "posedit length_change: return None on error"
        for tv in self.error_tests:
            v = self.hp.parse_hgvs_variant(tv)
            self.assertIsNone(v.posedit.length_change(on_error_raise=False))
