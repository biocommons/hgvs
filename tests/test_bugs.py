# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import pprint
import re
import sys
import unittest

from nose.plugins.attrib import attr

from hgvs.exceptions import HGVSError, HGVSDataNotAvailableError, HGVSParseError
import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.parser
import hgvs.transcriptmapper
import hgvs.variant
import hgvs.variantmapper



@attr(tags=["bugs"])
class Test_VariantMapper(unittest.TestCase):
    def setUp(self):
        self.hdp = hgvs.dataproviders.uta.connect()
        self.hm = hgvs.variantmapper.VariantMapper(self.hdp)
        self.hp = hgvs.parser.Parser()
        self.evm = hgvs.variantmapper.EasyVariantMapper(self.hdp,
                                                        replace_reference=True, primary_assembly='GRCh37',
                                                        alt_aln_method='splign')
        self.vn = hgvs.normalizer.Normalizer(self.hdp, shuffle_direction=3, cross_boundaries=True)


    def test_285_partial_palindrome_inversion(self):
        # https://bitbucket.org/biocommons/hgvs/issues/285/
        # Inversion mapping is not palindrome aware
        v = self.hp.parse_hgvs_variant("NM_000088.3:c.589_600inv")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NM_000088.3:c.590_599inv")


    def test_324_error_normalizing_simple_inversion(self):
        v = self.hp.parse_hgvs_variant("NM_000535.5:c.1673_1674inv")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NM_000535.5:c.1673_1674inv")  # no change

        vg = self.evm.c_to_g(v)
        v2 = self.evm.g_to_c(vg, tx_ac = v.ac)
        self.assertEqual(str(v), str(v2))  # no change after roundtrip


    def test_334_delins_normalization(self):
        # also tests 335 re: inv including sequence (e.g., NC_000009.11:g.36233991_36233992invCA)
        v = self.hp.parse_hgvs_variant("NC_000009.11:g.36233991_36233992delCAinsAC")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NC_000009.11:g.36233991_36233992delCAinsAC")

        v = self.hp.parse_hgvs_variant("NM_000535.5:c.1673_1674delCCinsGG")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NM_000535.5:c.1673_1674inv")
        
        v = self.hp.parse_hgvs_variant("NM_000535.5:c.1_3delAGCinsGCT")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NM_000535.5:c.1_3inv")


    def test_346_reject_partial_alignments(self):
        # hgvs-346: verify that alignment data covers full-length transcript
        with self.assertRaises(HGVSDataNotAvailableError):
            hgvs.transcriptmapper.TranscriptMapper(self.hdp,
                                                   tx_ac="NM_001290223.1",
                                                   alt_ac="NC_000010.10",
                                                   alt_aln_method="splign")
                                                   

    def test_386_reject_discontiguous_alignments(self):
        # hgvs-386: raise exception when alignments are discontiguous
        self.evm.c_to_g(self.hp.parse_hgvs_variant("NM_001637.3:c.1582G>A"))  # control
        with self.assertRaises(HGVSDataNotAvailableError):
            self.evm.c_to_g(self.hp.parse_hgvs_variant("NM_033517.1:c.1403A>C"))
