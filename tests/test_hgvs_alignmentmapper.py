# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import unittest

import pytest
from support import CACHE

import hgvs
import hgvs.dataproviders.uta
import hgvs.location
import hgvs.parser
from hgvs.alignmentmapper import AlignmentMapper
from hgvs.enums import Datum
from hgvs.exceptions import (HGVSDataNotAvailableError, HGVSError,
                             HGVSInvalidIntervalError)


@pytest.mark.quick
class Test_AlignmentMapper(unittest.TestCase):
    ref = "GRCh37.p10"

    @classmethod
    def setUpClass(cls):
        cls.hdp = hgvs.dataproviders.uta.connect(mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE)
        cls.parser = hgvs.parser.Parser()

    def test_alignmentmapper_failures(self):
        hgvs.global_config.mapping.strict_bounds = True
        with self.assertRaises(HGVSDataNotAvailableError):
            AlignmentMapper(self.hdp, tx_ac="bogus", alt_ac="NM_033089.6", alt_aln_method="splign")
        with self.assertRaises(HGVSDataNotAvailableError):
            AlignmentMapper(self.hdp, tx_ac="bogus", alt_ac="NM_033089.6", alt_aln_method="transcript")
        with self.assertRaises(HGVSDataNotAvailableError):
            AlignmentMapper(self.hdp, tx_ac="NM_033089.6", alt_ac="bogus", alt_aln_method="splign")
        with self.assertRaises(HGVSDataNotAvailableError):
            AlignmentMapper(self.hdp, tx_ac="NM_000051.3", alt_ac="NC_000011.9", alt_aln_method="bogus")
        with self.assertRaises(HGVSInvalidIntervalError):
            AlignmentMapper(self.hdp, "NM_000348.3", "NC_000002.11", "splign").n_to_c(
                self.parser.parse_n_interval("-1")
            )
        with self.assertRaises(HGVSInvalidIntervalError):
            AlignmentMapper(self.hdp, "NM_000348.3", "NC_000002.11", "splign").c_to_n(
                self.parser.parse_c_interval("99999")
            )

    def x_test_alignmentmapper_AlignmentMapper_LCE3C_uncertain(self):
        # ? is not yet supported
        """Use NM_178434.2 tests to test mapping with uncertain positions"""
        tx_ac = "NM_178434.2"
        alt_ac = "NC_000001.10"
        test_cases = [
            {
                "g": parser.parse_g_interval("(?_152573139)"),
                "n": parser.parse_n_interval("(?_2)"),
                "c": parser.parse_c_interval("(?_-69)"),
            },
            {
                "g": parser.parse_g_interval("(152573138_?)"),
                "n": parser.parse_n_interval("(1_?)"),
                "c": parser.parse_c_interval("(-70_?)"),
            },
        ]
        self.run_cases(tx_ac, alt_ac, test_cases)

    def test_alignmentmapper_AlignmentMapper_LCE3C(self):
        """NM_178434.2: LCE3C single exon, strand = +1, all coordinate input/output are in HGVS"""
        tx_ac = "NM_178434.2"
        alt_ac = "NC_000001.10"
        test_cases = [
            # 5'
            {
                "g": self.parser.parse_g_interval("152573138"),
                "n": self.parser.parse_n_interval("1"),
                "c": self.parser.parse_c_interval("-70"),
            },
            {
                "g": self.parser.parse_g_interval("152573140"),
                "n": self.parser.parse_n_interval("3"),
                "c": self.parser.parse_c_interval("-68"),
            },
            # cds
            {
                "g": self.parser.parse_g_interval("152573207"),
                "n": self.parser.parse_n_interval("70"),
                "c": self.parser.parse_c_interval("-1"),
            },
            {
                "g": self.parser.parse_g_interval("152573208"),
                "n": self.parser.parse_n_interval("71"),
                "c": self.parser.parse_c_interval("1"),
            },
            # 3'
            {
                "g": self.parser.parse_g_interval("152573492"),
                "n": self.parser.parse_n_interval("355"),
                "c": self.parser.parse_c_interval("285"),
            },
            {
                "g": self.parser.parse_g_interval("152573493"),
                "n": self.parser.parse_n_interval("356"),
                "c": self.parser.parse_c_interval("*1"),
            },
            {
                "g": self.parser.parse_g_interval("152573560"),
                "n": self.parser.parse_n_interval("423"),
                "c": self.parser.parse_c_interval("*68"),
            },
            {
                "g": self.parser.parse_g_interval("152573562"),
                "n": self.parser.parse_n_interval("425"),
                "c": self.parser.parse_c_interval("*70"),
            },
        ]
        self.run_cases(tx_ac, alt_ac, test_cases)

    def test_alignmentmapper_AlignmentMapper_HIST3H2A(self):
        """NM_033445.2: LCE3C single exon, strand = -1, all coordinate input/output are in HGVS"""
        tx_ac = "NM_033445.2"
        alt_ac = "NC_000001.10"
        test_cases = [
            # 3'
            {
                "g": self.parser.parse_g_interval("228645560"),
                "n": self.parser.parse_n_interval("1"),
                "c": self.parser.parse_c_interval("-42"),
            },
            {
                "g": self.parser.parse_g_interval("228645558"),
                "n": self.parser.parse_n_interval("3"),
                "c": self.parser.parse_c_interval("-40"),
            },
            # cds
            {
                "g": self.parser.parse_g_interval("228645519"),
                "n": self.parser.parse_n_interval("42"),
                "c": self.parser.parse_c_interval("-1"),
            },
            {
                "g": self.parser.parse_g_interval("228645518"),
                "n": self.parser.parse_n_interval("43"),
                "c": self.parser.parse_c_interval("1"),
            },
            # 5'
            {
                "g": self.parser.parse_g_interval("228645126"),
                "n": self.parser.parse_n_interval("435"),
                "c": self.parser.parse_c_interval("393"),
            },
            {
                "g": self.parser.parse_g_interval("228645125"),
                "n": self.parser.parse_n_interval("436"),
                "c": self.parser.parse_c_interval("*1"),
            },
            {
                "g": self.parser.parse_g_interval("228645124"),
                "n": self.parser.parse_n_interval("437"),
                "c": self.parser.parse_c_interval("*2"),
            },
            {
                "g": self.parser.parse_g_interval("228645065"),
                "n": self.parser.parse_n_interval("496"),
                "c": self.parser.parse_c_interval("*61"),
            },
        ]
        self.run_cases(tx_ac, alt_ac, test_cases)

    def test_alignmentmapper_AlignmentMapper_LCE2B(self):
        """NM_014357.4: LCE2B, two exons, strand = +1, all coordinate input/output are in HGVS"""
        tx_ac = "NM_014357.4"
        alt_ac = "NC_000001.10"
        test_cases = [
            # 5'
            {
                "g": self.parser.parse_g_interval("152658599"),
                "n": self.parser.parse_n_interval("1"),
                "c": self.parser.parse_c_interval("-54"),
            },
            {
                "g": self.parser.parse_g_interval("152658601"),
                "n": self.parser.parse_n_interval("3"),
                "c": self.parser.parse_c_interval("-52"),
            },
            # cds
            {
                "g": self.parser.parse_g_interval("152659319"),
                "n": self.parser.parse_n_interval("54"),
                "c": self.parser.parse_c_interval("-1"),
            },
            {
                "g": self.parser.parse_g_interval("152659320"),
                "n": self.parser.parse_n_interval("55"),
                "c": self.parser.parse_c_interval("1"),
            },
            # around end of exon 1
            {
                "g": self.parser.parse_g_interval("152658632"),
                "n": self.parser.parse_n_interval("34"),
                "c": self.parser.parse_c_interval("-21"),
            },
            {
                "g": self.parser.parse_g_interval("152658633"),
                "n": self.parser.parse_n_interval("34+1"),
                "c": self.parser.parse_c_interval("-21+1"),
            },
            # span
            {
                "g": self.parser.parse_g_interval("152658633_152659299"),
                "n": self.parser.parse_n_interval("34+1_35-1"),
                "c": self.parser.parse_c_interval("-21+1_-20-1"),
            },
            # around beginning of exon 2
            {
                "g": self.parser.parse_g_interval("152659300"),
                "n": self.parser.parse_n_interval("35"),
                "c": self.parser.parse_c_interval("-20"),
            },
            {
                "g": self.parser.parse_g_interval("152659299"),
                "n": self.parser.parse_n_interval("35-1"),
                "c": self.parser.parse_c_interval("-20-1"),
            },
            # around end of exon 2
            {
                "g": self.parser.parse_g_interval("152659652"),
                "n": self.parser.parse_n_interval("387"),
                "c": self.parser.parse_c_interval("333"),
            },
            {
                "g": self.parser.parse_g_interval("152659653"),
                "n": self.parser.parse_n_interval("388"),
                "c": self.parser.parse_c_interval("*1"),
            },
            # span
            {
                "g": self.parser.parse_g_interval("152659651_152659654"),
                "n": self.parser.parse_n_interval("386_389"),
                "c": self.parser.parse_c_interval("332_*2"),
            },
            # 3'
            {
                "g": self.parser.parse_g_interval("152659877"),
                "n": self.parser.parse_n_interval("612"),
                "c": self.parser.parse_c_interval("*225"),
            },
        ]
        self.run_cases(tx_ac, alt_ac, test_cases)

    def test_alignmentmapper_AlignmentMapper_PTH2(self):
        """NM_178449.3: PTH2, two exons, strand = -1, all coordinate input/output are in HGVS"""
        tx_ac = "NM_178449.3"
        alt_ac = "NC_000019.9"
        test_cases = [
            # 3'
            {
                "g": self.parser.parse_g_interval("49926698"),
                "n": self.parser.parse_n_interval("1"),
                "c": self.parser.parse_c_interval("-102"),
            },
            # cds
            {
                "g": self.parser.parse_g_interval("49926597"),
                "n": self.parser.parse_n_interval("102"),
                "c": self.parser.parse_c_interval("-1"),
            },
            {
                "g": self.parser.parse_g_interval("49926596"),
                "n": self.parser.parse_n_interval("103"),
                "c": self.parser.parse_c_interval("1"),
            },
            # around end of exon 1
            {
                "g": self.parser.parse_g_interval("49926469"),
                "n": self.parser.parse_n_interval("230"),
                "c": self.parser.parse_c_interval("128"),
            },
            {
                "g": self.parser.parse_g_interval("49926468"),
                "n": self.parser.parse_n_interval("230+1"),
                "c": self.parser.parse_c_interval("128+1"),
            },
            # span
            {
                "g": self.parser.parse_g_interval("49925901_49926467"),
                "n": self.parser.parse_n_interval("230+2_231-2"),
                "c": self.parser.parse_c_interval("128+2_129-2"),
            },
            # around beginning of exon 2
            {
                "g": self.parser.parse_g_interval("49925900"),
                "n": self.parser.parse_n_interval("231-1"),
                "c": self.parser.parse_c_interval("129-1"),
            },
            {
                "g": self.parser.parse_g_interval("49925899"),
                "n": self.parser.parse_n_interval("231"),
                "c": self.parser.parse_c_interval("129"),
            },
            # around end of exon 2
            {
                "g": self.parser.parse_g_interval("49925725"),
                "n": self.parser.parse_n_interval("405"),
                "c": self.parser.parse_c_interval("303"),
            },
            {
                "g": self.parser.parse_g_interval("49925724"),
                "n": self.parser.parse_n_interval("406"),
                "c": self.parser.parse_c_interval("*1"),
            },
            {
                "g": self.parser.parse_g_interval("49925671"),
                "n": self.parser.parse_n_interval("459"),
                "c": self.parser.parse_c_interval("*54"),
            },
        ]
        self.run_cases(tx_ac, alt_ac, test_cases)

    def run_cases(self, tx_ac, alt_ac, test_cases):
        am = AlignmentMapper(self.hdp, tx_ac, alt_ac, alt_aln_method="splign")
        for test_case in test_cases:
            assert test_case["c"] == am.g_to_c(test_case["g"]), f"{tx_ac}~{alt_ac} {test_case['c']} am.g_to_c failed"
            assert test_case["c"] == am.n_to_c(test_case["n"]), f"{tx_ac}~{alt_ac} {test_case['c']} am.n_to_c failed"
            assert test_case["g"] == am.c_to_g(test_case["c"]), f"{tx_ac}~{alt_ac} {test_case['g']} am.c_to_g failed"
            assert test_case["g"] == am.n_to_g(test_case["n"]), f"{tx_ac}~{alt_ac} {test_case['g']} am.n_to_g failed"
            assert test_case["n"] == am.c_to_n(test_case["c"]), f"{tx_ac}~{alt_ac} {test_case['n']} am.c_to_n failed"
            assert test_case["n"] == am.g_to_n(test_case["g"]), f"{tx_ac}~{alt_ac} {test_case['n']} am.g_to_n failed"


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
