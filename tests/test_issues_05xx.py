# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import pprint
import re
import sys
import os

import unittest

import pytest

from hgvs.exceptions import HGVSError, HGVSDataNotAvailableError, HGVSParseError, HGVSInvalidVariantError
from hgvs.enums import Datum
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.parser
import hgvs.sequencevariant
import hgvs.validator
import hgvs.variantmapper
from support import CACHE


@pytest.mark.issues
class Test_Issues(unittest.TestCase):
    def setUp(self):
        self.hdp = hgvs.dataproviders.uta.connect(
            mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE)
        self.vm = hgvs.variantmapper.VariantMapper(self.hdp, replace_reference=False)
        self.vm_rr = hgvs.variantmapper.VariantMapper(self.hdp, replace_reference=True)
        self.hp = hgvs.parser.Parser()
        self.hn = hgvs.normalizer.Normalizer(self.hdp)
        self.hv = hgvs.validator.IntrinsicValidator()
        self.am37 = hgvs.assemblymapper.AssemblyMapper(
            self.hdp, replace_reference=True, assembly_name='GRCh37', alt_aln_method='splign')
        self.am38 = hgvs.assemblymapper.AssemblyMapper(
            self.hdp, replace_reference=True, assembly_name='GRCh38', alt_aln_method='splign')
        self.vn = hgvs.normalizer.Normalizer(self.hdp, shuffle_direction=3, cross_boundaries=True)

    def test_525(self):
        """https://github.com/biocommons/hgvs/issues/525"""

        # simple test case
        hgvs = "NM_001637.3:c.3_4insTAG"    # insert stop in phase at AA 2
        var_c = self.hp.parse_hgvs_variant(hgvs)
        var_p = self.am38.c_to_p(var_c)
        assert str(var_p) == "NP_001628.1:p.(Gln2Ter)"
