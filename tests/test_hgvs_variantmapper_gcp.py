# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import pprint
import re
import sys
import unittest

import unicodecsv as csv

from nose.plugins.attrib import attr

from hgvs.exceptions import HGVSError
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variant
import hgvs.variantmapper


def gxp_file_reader(fn):
    rdr = csv.DictReader(open(fn, 'r'), delimiter=str('\t'))
    for rec in rdr:
        if rec['id'].startswith('#'):
            continue
        yield rec


@attr(tags=["mapping"])
class Test_VariantMapper(unittest.TestCase):
    def setUp(self):
        self.hdp = hgvs.dataproviders.uta.connect()
        self.hm = hgvs.variantmapper.VariantMapper(self.hdp)
        self.hp = hgvs.parser.Parser()

    # ZCCHC3 -- one exon, + strand
    # reece@[local]/uta_dev=> select hgnc,alt_strand,n_exons,tx_ac,alt_ac,s_cigars,cds_start_i,cds_end_i from bermuda.bermuda_data_mv where tx_ac = 'NM_033089.6';
    # ┌────────┬────────────┬─────────┬─────────────┬──────────────┬─────────────┬─────────────┬───────────┐
    # │  hgnc  │ alt_strand │ n_exons │    tx_ac    │    alt_ac    │  s_cigars   │ cds_start_i │ cds_end_i │
    # ├────────┼────────────┼─────────┼─────────────┼──────────────┼─────────────┼─────────────┼───────────┤
    # │ ZCCHC3 │          1 │       1 │ NM_033089.6 │ NC_000020.10 │ 484=3I2275= │          24 │      1236 │
    # └────────┴────────────┴─────────┴─────────────┴──────────────┴─────────────┴─────────────┴───────────┘
    def test_ZCCHC3_dbSNP(self):
        for rec in gxp_file_reader('tests/data/gcp/ZCCHC3-dbSNP.tsv'):
            self._test_gxp_mapping(rec)

    # ORAI1 -- two exons, + strand
    # reece@[local]/uta_dev=> select hgnc,alt_strand,n_exons,tx_ac,alt_ac,s_cigars,cds_start_i,cds_end_i from bermuda.bermuda_data_mv where tx_ac = 'NM_032790.3';
    # ┌───────┬────────────┬─────────┬─────────────┬──────────────┬──────────────────┬─────────────┬───────────┐
    # │ hgnc  │ alt_strand │ n_exons │    tx_ac    │    alt_ac    │     s_cigars     │ cds_start_i │ cds_end_i │
    # ├───────┼────────────┼─────────┼─────────────┼──────────────┼──────────────────┼─────────────┼───────────┤
    # │ ORAI1 │          1 │       2 │ NM_032790.3 │ NC_000012.11 │ 319=6I177=;1000= │         193 │      1099 │
    # └───────┴────────────┴─────────┴─────────────┴──────────────┴──────────────────┴─────────────┴───────────┘
    def test_ORAI1_dbSNP(self):
        for rec in gxp_file_reader('tests/data/gcp/ORAI1-dbSNP.tsv'):
            self._test_gxp_mapping(rec)

    # FOLR3 -- multiple exons, + strand
    # reece@[local]/uta_dev=> select hgnc,alt_strand,n_exons,tx_ac,alt_ac,s_cigars,cds_start_i,cds_end_i from bermuda.bermuda_data_mv where tx_ac = 'NM_000804.2';
    # ┌───────┬────────────┬─────────┬─────────────┬─────────────┬──────────────────────────────┬─────────────┬───────────┐
    # │ hgnc  │ alt_strand │ n_exons │    tx_ac    │   alt_ac    │           s_cigars           │ cds_start_i │ cds_end_i │
    # ├───────┼────────────┼─────────┼─────────────┼─────────────┼──────────────────────────────┼─────────────┼───────────┤
    # │ FOLR3 │          1 │       5 │ NM_000804.2 │ NC_000011.9 │ 44=;174=;150=2D37=;136=;304= │          50 │       788 │
    # └───────┴────────────┴─────────┴─────────────┴─────────────┴──────────────────────────────┴─────────────┴───────────┘
    def test_FOLR3_dbSNP(self):
        # TODO: CORE-158: g-to-c mapped insertions have incorrect interval bounds
        for rec in gxp_file_reader('tests/data/gcp/FOLR3-dbSNP.tsv'):
            self._test_gxp_mapping(rec)

    # ADRA2B -- one exon, - strand
    # reece@[local]/uta_dev=> select hgnc,alt_strand,n_exons,tx_ac,alt_ac,s_cigars,cds_start_i,cds_end_i from bermuda.bermuda_data_mv where tx_ac = 'NM_000682.5';
    # ┌────────┬────────────┬─────────┬─────────────┬──────────────┬─────────────┬─────────────┬───────────┐
    # │  hgnc  │ alt_strand │ n_exons │    tx_ac    │    alt_ac    │  s_cigars   │ cds_start_i │ cds_end_i │
    # ├────────┼────────────┼─────────┼─────────────┼──────────────┼─────────────┼─────────────┼───────────┤
    # │ ADRA2B │         -1 │       1 │ NM_000682.5 │ NC_000002.11 │ 891=9D2375= │           0 │      1353 │
    # └────────┴────────────┴─────────┴─────────────┴──────────────┴─────────────┴─────────────┴───────────┘
    def test_ADRA2B_dbSNP(self):
        for rec in gxp_file_reader('tests/data/gcp/ADRA2B-dbSNP.tsv'):
            self._test_gxp_mapping(rec)

    # JRK -- multiple exons, - strand
    # reece@[local]/uta_dev=> select hgnc,alt_strand,n_exons,tx_ac,alt_ac,s_cigars,cds_start_i,cds_end_i from bermuda.bermuda_data_mv where tx_ac = 'NM_001077527.1';
    # ┌──────┬────────────┬─────────┬────────────────┬──────────────┬───────────────────────┬─────────────┬───────────┐
    # │ hgnc │ alt_strand │ n_exons │     tx_ac      │    alt_ac    │       s_cigars        │ cds_start_i │ cds_end_i │
    # ├──────┼────────────┼─────────┼────────────────┼──────────────┼───────────────────────┼─────────────┼───────────┤
    # │ JRK  │         -1 │       3 │ NM_001077527.1 │ NC_000008.10 │ 52=;1844=2I199=;1483= │         514 │      2185 │
    # └──────┴────────────┴─────────┴────────────────┴──────────────┴───────────────────────┴─────────────┴───────────┘
    def test_JRK_dbSNP(self):
        # TODO: CORE-157: del26 on -1 strands gets reverse complemented as del62
        for rec in gxp_file_reader('tests/data/gcp/JRK-dbSNP.tsv'):
            self._test_gxp_mapping(rec)

    def test_NEFL_dbSNP(self):
        for rec in gxp_file_reader('tests/data/gcp/NEFL-dbSNP.tsv'):
            self._test_gxp_mapping(rec)

    def test_DNAH11_hgmd(self):
        for rec in gxp_file_reader('tests/data/gcp/DNAH11-HGMD.tsv'):
            self._test_gxp_mapping(rec)

    def test_DNAH11_dbSNP_NM_003777(self):
        for rec in gxp_file_reader('tests/data/gcp/DNAH11-dbSNP-NM_003777.tsv'):
            self._test_gxp_mapping(rec)

    def test_DNAH11_dbSNP_NM_001277115(self):
        for rec in gxp_file_reader('tests/data/gcp/DNAH11-dbSNP-NM_001277115.tsv'):
            self._test_gxp_mapping(rec)

    @attr(tags=["regression"])
    def test_regression(self):
        for rec in gxp_file_reader('tests/data/gcp/regression.tsv'):
            self._test_gxp_mapping(rec)

    @attr(tags=["extra"])
    def test_DNAH11_dbSNP_full(self):
        for rec in gxp_file_reader('tests/data/gcp/DNAH11-dbSNP.tsv'):
            self._test_gxp_mapping(rec)

    def test_real(self):
        for rec in gxp_file_reader('tests/data/gcp/real.tsv'):
            self._test_gxp_mapping(rec)

    def test_noncoding(self):
        for rec in gxp_file_reader('tests/data/gcp/noncoding.tsv'):
            self._test_gxp_mapping(rec)

    def _test_gxp_mapping(self, rec):
        """given one record (row) of g, c/n/r, and p (optional) test variants, map
        g->c/n/r, c/n/r->g, and c->p and verify equivalence

        """

        def _rm_del_seq(vs):
            return re.sub(vs, 'del\w+ins', 'delins')

        var_g = self.hp.parse_hgvs_variant(rec['HGVSg'])
        var_x = self.hp.parse_hgvs_variant(rec['HGVSc'])
        var_p = self.hp.parse_hgvs_variant(rec['HGVSp']) if rec['HGVSp'] is not None and rec['HGVSp'] != '' else None

        # g -> x
        if var_x.type == 'c':
            var_x_test = self.hm.g_to_c(var_g, var_x.ac)
        elif var_x.type == 'n':
            var_x_test = self.hm.g_to_n(var_g, var_x.ac)
        self.assertEquals(_rm_del_seq(str(var_x)), _rm_del_seq(str(var_x_test)),
                          msg="%s != %s (%s; HGVSg=%s)" % (str(var_x_test), str(var_x), rec['id'], rec['HGVSg']))

        # c,n -> g
        if var_x.type == 'c':
            var_g_test = self.hm.c_to_g(var_x, var_g.ac)
        elif var_x.type == 'n':
            var_g_test = self.hm.n_to_g(var_x, var_g.ac)
        self.assertEquals(_rm_del_seq(str(var_g)), _rm_del_seq(str(var_g_test)),
                          msg="%s != %s (%s; HGVSc=%s)" % (str(var_g_test), str(var_g), rec['id'], rec['HGVSc']))

        if var_p is not None:
            # c -> p
            hgvs_p_exp = str(var_p)
            var_p_test = self.hm.c_to_p(var_x, var_p.ac)

            if not var_p.posedit.uncertain:
                # if expected value isn't uncertain, strip uncertain from test
                var_p_test.posedit.uncertain = False

            hgvs_p_test = str(var_p_test)

            if re.search('Ter$', hgvs_p_exp):
                # if expected value doesn't have a count, strip it from the test
                hgvs_p_test = re.sub('Ter\d+$', 'Ter', hgvs_p_test)

            self.assertEquals(hgvs_p_exp, hgvs_p_test, msg="%s != %s (%s)" % (hgvs_p_exp, hgvs_p_test, rec['id']))


if __name__ == '__main__':
    unittest.main()

## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
## 
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
## 
##     http://www.apache.org/licenses/LICENSE-2.0
## 
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
