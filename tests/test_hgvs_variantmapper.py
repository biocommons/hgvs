# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import os

import unittest

import pytest

from hgvs.exceptions import HGVSError, HGVSInvalidVariantError
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper
from support import CACHE


def test_add_gene_symbol(am38, parser):
    ags = am38.add_gene_symbol
    var_g = parser.parse("NC_000007.13:g.21940852_21940908del")

    am38.add_gene_symbol = True
    var_t = am38.g_to_t(var_g, "NM_003777.3")
    assert "NM_003777.3(DNAH11):c.13552_*36del" == str(var_t)

    am38.add_gene_symbol = False
    var_t = am38.g_to_t(var_g, "NM_003777.3")
    assert "NM_003777.3:c.13552_*36del" == str(var_t)

    am38.add_gene_symbol = ags


@pytest.mark.quick
class Test_VariantMapper_Exceptions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hdp = hgvs.dataproviders.uta.connect(
            mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE)
        cls.vm = hgvs.variantmapper.VariantMapper(cls.hdp)
        cls.hp = hgvs.parser.Parser()

    def test_gcrp_invalid_input_type(self):
        hgvs_g = "NC_000007.13:g.36561662C>T"
        hgvs_c = "NM_001637.3:c.1582G>A"

        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_c = self.hp.parse_hgvs_variant(hgvs_c)

        cases = {
            "gc": (self.vm.g_to_c, (var_c, "NM_001637.3")),
            "gr": (self.vm.g_to_n, (var_c, "NM_001637.3")),
            "gt": (self.vm.g_to_t, (var_c, "NM_001637.3")),
            "rg": (self.vm.n_to_g, (var_c, "NM_001637.3")),
            "cg": (self.vm.c_to_g, (var_g, "NM_001637.3")),
            "tg": (self.vm.t_to_g, (var_g, "NM_001637.3")),
            "cr": (self.vm.c_to_n, (var_g, )),
            "rc": (self.vm.n_to_c, (var_g, )),
            "cp": (self.vm.c_to_p, (var_g, None)),
        }

        failures = []
        for key in cases:
            try:
                func, args = cases[key]
                var_result = func(*args)
                failures.append(key)
            except hgvs.exceptions.HGVSInvalidVariantError:
                pass

        self.assertFalse(failures, "conversions not failing: {}".format(failures))

    def test_gc_invalid_input_nm_accession(self):
        hgvs_g = "NC_000007.13:g.36561662C>T"
        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        with self.assertRaises(hgvs.exceptions.HGVSError):
            var_p = self.vm.c_to_p(var_g, "NM_999999.1")

    def test_undefined_cds(self):
        """Raise exception when requesting mapping to/from c. with non-coding transcript"""
        hgvs_n = "NR_111984.1:n.44G>A"    # legit
        hgvs_c = "NR_111984.1:c.44G>A"    # bogus: c. with non-coding tx accession
        var_n = self.hp.parse_hgvs_variant(hgvs_n)
        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        tx_ac = var_n.ac

        with self.assertRaises(hgvs.exceptions.HGVSUsageError):
            var_c = self.vm.n_to_c(var_n)    # n_to_c: transcript is non-coding

        with self.assertRaises(hgvs.exceptions.HGVSUsageError):
            var_c = self.vm.c_to_n(var_c)    # c_to_n: var_c is bogus

    def test_map_var_of_unsupported_validation(self):
        hgvs_c = "NM_003777.3:c.13552_*36del57"
        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_g = self.vm.c_to_g(var_c, "NC_000007.13")
        self.assertEqual(str(var_g), "NC_000007.13:g.21940852_21940908del")

    def test_map_to_unknown_p_effect(self):
        hgvs_c = "NM_020975.4:c.625+9C>T"
        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.vm.c_to_p(var_c)
        self.assertEqual(str(var_p), "NP_066124.1:p.?")

    def test_map_of_c_out_of_cds_bound(self):
        hgvs_c = "NM_145901.2:c.343T>C"
        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        with self.assertRaises(HGVSInvalidVariantError):
            self.vm.c_to_p(var_c)

    def test_map_of_dup_at_cds_end(self):
        hgvs_c = "NM_001051.2:c.1257dupG"
        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.vm.c_to_p(var_c)
        self.assertEqual(str(var_p), "NP_001042.1:p.(=)")

    def test_map_of_c_out_of_reference_bound(self):
        hgvs_c = "NM_000249.3:c.-73960_*46597del"
        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        with pytest.raises(HGVSError, match='coordinate is outside the bounds'):
            self.vm.c_to_p(var_c)


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
