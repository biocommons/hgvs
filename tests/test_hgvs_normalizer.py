# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import os

import unittest

import pytest

from hgvs.exceptions import HGVSError, HGVSUnsupportedOperationError, HGVSInvalidVariantError
import hgvs.dataproviders.uta
import hgvs.variantmapper
import hgvs.parser
import hgvs.normalizer
from support import CACHE

hdp = hgvs.dataproviders.uta.connect(mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE)


@pytest.mark.normalization
class Test_HGVSNormalizer(unittest.TestCase):
    """Tests for normalizer"""

    @classmethod
    def setUpClass(cls):
        cls.hp = hgvs.parser.Parser()
        cls.norm = hgvs.normalizer.Normalizer(hdp, shuffle_direction=3, cross_boundaries=True)
        cls.norm5 = hgvs.normalizer.Normalizer(hdp, shuffle_direction=5, cross_boundaries=True)
        cls.normc = hgvs.normalizer.Normalizer(hdp, shuffle_direction=3, cross_boundaries=False)
        cls.norm5c = hgvs.normalizer.Normalizer(hdp, shuffle_direction=5, cross_boundaries=False)

    def test_c_normalizer(self):
        """Test normalizer for variant type c."""
        #3' shuffling
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_000088.3:c.589_600inv"))),
            "NM_000088.3:c.590_599inv")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.31del"))),
            "NM_001166478.1:c.35del")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.35_36insT"))),
            "NM_001166478.1:c.35dup")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.36_37insTC"))),
            "NM_001166478.1:c.36_37dup")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.35_36dup"))),
            "NM_001166478.1:c.36_37dup")
        self.assertEqual(
            str(
                self.norm.normalize(
                    self.hp.parse_hgvs_variant("NM_001166478.1:c.2_7delinsTTTAGA"))),
            "NM_001166478.1:c.3_4delinsTT")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.30_31insT"))),
            "NM_001166478.1:c.35dup")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.59delG"))),
            "NM_001166478.1:c.61del")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.36_37insTCTCTC"))),
            "NM_001166478.1:c.37_38insCTCTCT")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_000051.3:c.14_15insT"))),
            "NM_000051.3:c.15dup")

        #5' shuffling
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_000088.3:c.589_600inv"))),
            "NM_000088.3:c.590_599inv")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.34del"))),
            "NM_001166478.1:c.31del")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.35_36insT"))),
            "NM_001166478.1:c.31dup")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.36_37insTC"))),
            "NM_001166478.1:c.35_36dup")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.35_36dup"))),
            "NM_001166478.1:c.35_36dup")
        self.assertEqual(
            str(
                self.norm5.normalize(
                    self.hp.parse_hgvs_variant("NM_001166478.1:c.2_7delinsTTTAGA"))),
            "NM_001166478.1:c.3_4delinsTT")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.30_31insT"))),
            "NM_001166478.1:c.31dup")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.61delG"))),
            "NM_001166478.1:c.59del")
        self.assertEqual(
            str(
                self.norm5.normalize(
                    self.hp.parse_hgvs_variant("NM_001166478.1:c.36_37insTCTCTC"))),
            "NM_001166478.1:c.34_35insTCTCTC")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_000051.3:c.14_15insT"))),
            "NM_000051.3:c.14dup")

        #Around exon-intron boundary
        self.assertEqual(
            str(self.normc.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.59delG"))),
            "NM_001166478.1:c.60del")
        self.assertEqual(
            str(self.norm5c.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.61delG"))),
            "NM_001166478.1:c.61del")
        self.assertEqual(
            str(self.norm5c.normalize(self.hp.parse_hgvs_variant("NM_001110792.1:c.1030_1035del"))),
            "NM_001110792.1:c.1029_1034del")
        with self.assertRaises(HGVSUnsupportedOperationError):
            self.normc.normalize(self.hp.parse_hgvs_variant("NM_001166478.1:c.59_61del"))

        #UTR variants
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_000051.3:c.-5_-4insA"))),
            "NM_000051.3:c.-3dup")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_000051.3:c.-4_-3insAC"))),
            "NM_000051.3:c.-3_-2dup")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_000051.3:c.-2_-1insCA"))),
            "NM_000051.3:c.-1_1dup")

        self.assertEqual(
            str(self.normc.normalize(self.hp.parse_hgvs_variant("NM_000051.3:c.-2_-1insCA"))),
            "NM_000051.3:c.-1_1insAC")

        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_000051.3:c.-4_-3insA"))),
            "NM_000051.3:c.-4dup")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_000051.3:c.1_2insCA"))),
            "NM_000051.3:c.-1_1dup")

        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_000051.3:c.*2_*3insT"))),
            "NM_000051.3:c.*4dup")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_000051.3:c.9170_9171insAT"))),
            "NM_000051.3:c.9171_*1dup")

        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_000051.3:c.*4_*5insT"))),
            "NM_000051.3:c.*3dup")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_000051.3:c.9171_*1insA"))),
            "NM_000051.3:c.9171dup")

        with self.assertRaises(HGVSInvalidVariantError):
            self.norm.normalize(self.hp.parse_hgvs_variant("NM_000059.3:c.7790delAAG"))

    def test_g_normalizer(self):
        """Test normalizer for variant type g."""
        #3' shuffling
        self.assertEqual(
            str(
                self.norm.normalize(
                    self.hp.parse_hgvs_variant("NC_000006.11:g.49917122_49917123insA"))),
            "NC_000006.11:g.49917127dup")
        self.assertEqual(
            str(
                self.norm.normalize(
                    self.hp.parse_hgvs_variant("NC_000006.11:g.49917121_49917122insGA"))),
            "NC_000006.11:g.49917122_49917123dup")
        self.assertEqual(
            str(
                self.norm.normalize(
                    self.hp.parse_hgvs_variant("NC_000006.11:g.49917122_49917123dup"))),
            "NC_000006.11:g.49917122_49917123dup")
        self.assertEqual(
            str(
                self.norm.normalize(
                    self.hp.parse_hgvs_variant("NC_000006.11:g.49917122_49917123dupGA"))),
            "NC_000006.11:g.49917122_49917123dup")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NC_000006.11:g.49917098delC"))),
            "NC_000006.11:g.49917099del")
        self.assertEqual(
            str(
                self.norm.normalize(
                    self.hp.parse_hgvs_variant("NC_000006.11:g.49917151_49917156delinsTCTAAA"))),
            "NC_000006.11:g.49917154_49917155delinsAA")
        self.assertEqual(
            str(
                self.norm.normalize(
                    self.hp.parse_hgvs_variant("NC_000001.10:g.1647893delinsCTTTCTT"))),
            "NC_000001.10:g.1647895_1647900dup")
        self.assertEqual(
            str(
                self.norm.normalize(
                    self.hp.parse_hgvs_variant(
                        "NC_000003.12:g.46709584_46709610del27insAAGAAGAAGAAGAAGAAGAAGAAGAAG"))),
            "NC_000003.12:g.46709584_46709610=")

        #5' shuffling
        self.assertEqual(
            str(
                self.norm5.normalize(
                    self.hp.parse_hgvs_variant("NC_000006.11:g.49917122_49917123insA"))),
            "NC_000006.11:g.49917123dup")
        self.assertEqual(
            str(
                self.norm5.normalize(
                    self.hp.parse_hgvs_variant("NC_000006.11:g.49917121_49917122insGA"))),
            "NC_000006.11:g.49917121_49917122dup")
        self.assertEqual(
            str(
                self.norm5.normalize(
                    self.hp.parse_hgvs_variant("NC_000006.11:g.49917122_49917123dup"))),
            "NC_000006.11:g.49917121_49917122dup")
        self.assertEqual(
            str(
                self.norm5.normalize(
                    self.hp.parse_hgvs_variant("NC_000006.11:g.49917122_49917123dupGA"))),
            "NC_000006.11:g.49917121_49917122dup")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NC_000006.11:g.49917099delC"))),
            "NC_000006.11:g.49917098del")
        self.assertEqual(
            str(
                self.norm5.normalize(
                    self.hp.parse_hgvs_variant("NC_000006.11:g.49917151_49917156delinsTCTAAA"))),
            "NC_000006.11:g.49917154_49917155delinsAA")
        self.assertEqual(
            str(
                self.norm5.normalize(
                    self.hp.parse_hgvs_variant(
                        "NC_000003.12:g.46709584_46709610del27insAAGAAGAAGAAGAAGAAGAAGAAGAAG"))),
            "NC_000003.12:g.46709584_46709610=")

        self.assertEqual(
            str(
                self.norm.normalize(
                    self.hp.parse_hgvs_variant("NC_000009.11:g.36233991_36233992delCAinsTG"))),
            "NC_000009.11:g.36233991_36233992inv")
        with self.assertRaises(HGVSInvalidVariantError):
            self.norm.normalize(
                self.hp.parse_hgvs_variant("NG_032871.1:g.32476_53457delinsAATTAAGGTATA"))

        with self.assertRaises(HGVSInvalidVariantError):
            self.norm.normalize(
                self.hp.parse_hgvs_variant("NG_032871.1:g.32476_53457delinsAATTAAGGTATA"))

    def test_norm_var_near_bound(self):
        """Test normalizing variants near the end or start of transcript"""
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.935T>C"))),
            "NM_001001656.1:c.935T>C")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.945G>C"))),
            "NM_001001656.1:c.945G>C")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.945dup"))),
            "NM_001001656.1:c.945dup")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.935_945del"))),
            "NM_001001656.1:c.935_945del")
        with self.assertRaises(HGVSError):
            self.norm.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.946G>C"))
        with self.assertRaises(HGVSError):
            self.norm.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.946dup"))
        with self.assertRaises(HGVSError):
            self.norm.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.935_946del"))

        self.assertEqual(
            str(self.normc.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.935T>C"))),
            "NM_001001656.1:c.935T>C")
        self.assertEqual(
            str(self.normc.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.945G>C"))),
            "NM_001001656.1:c.945G>C")
        self.assertEqual(
            str(self.normc.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.945dup"))),
            "NM_001001656.1:c.945dup")
        self.assertEqual(
            str(self.normc.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.935_945del"))),
            "NM_001001656.1:c.935_945del")
        with self.assertRaises(HGVSError):
            self.normc.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.946G>C"))
        with self.assertRaises(HGVSError):
            self.normc.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.946dup"))
        with self.assertRaises(HGVSError):
            self.normc.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.935_946del"))

        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.1A>T"))),
            "NM_001001656.1:c.1A>T")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.1del"))),
            "NM_001001656.1:c.1del")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.1dup"))),
            "NM_001001656.1:c.1dup")

        self.assertEqual(
            str(self.norm5c.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.1A>T"))),
            "NM_001001656.1:c.1A>T")
        self.assertEqual(
            str(self.norm5c.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.1del"))),
            "NM_001001656.1:c.1del")
        self.assertEqual(
            str(self.norm5c.normalize(self.hp.parse_hgvs_variant("NM_001001656.1:c.1dup"))),
            "NM_001001656.1:c.1dup")

        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_212556.2:c.1_2insCA"))),
            "NM_212556.2:c.1_2insCA")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_212556.2:c.2_3insCAT"))),
            "NM_212556.2:c.2_3insCAT")
        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_212556.2:c.1delinsCA"))),
            "NM_212556.2:c.1delinsCA")

        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_212556.2:c.1_2insCA"))),
            "NM_212556.2:c.1delinsACA")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_212556.2:c.2_3insCAT"))),
            "NM_212556.2:c.1delinsATCA")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_212556.2:c.1delinsCA"))),
            "NM_212556.2:c.1delinsCA")
        self.assertEqual(
            str(self.norm5.normalize(self.hp.parse_hgvs_variant("NM_212556.2:c.1delinsAA"))),
            "NM_212556.2:c.1dup")

        self.assertEqual(
            str(self.norm.normalize(self.hp.parse_hgvs_variant("NM_212556.2:c.1400_1401insAC"))),
            "NM_212556.2:c.1401delinsACA")
        self.assertEqual(
            str(self.normc.normalize(self.hp.parse_hgvs_variant("NM_212556.2:c.1400_1401insAC"))),
            "NM_212556.2:c.1401delinsACA")
        self.assertEqual(
            str(self.normc.normalize(self.hp.parse_hgvs_variant("NM_212556.2:c.1401delinsAA"))),
            "NM_212556.2:c.1401dup")


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
