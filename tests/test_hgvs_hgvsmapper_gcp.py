## These tests are currently disabled because a large number of them still fail
# to execute, try:
# PYTHONPATH=. python tests/test_hgvs_hgvsmapper_gcp.py

import csv
import unittest

import bdi.sources.uta0_sqlite

import hgvs.hgvsmapper as hgvsmapper
import hgvs.parser


def gcp_file_reader(fn):
    rdr = csv.DictReader(open(fn, 'r'), delimiter='\t')
    for rec in rdr:
        if rec['id'].startswith('#'):
            continue
        yield rec


class Test_HGVSMapper(unittest.TestCase):
    def setUp(self):
        self.bdi = bdi.sources.uta0_sqlite.UTA0('/tmp/uta-0.0.4.db')
        self.hm = hgvs.hgvsmapper.HGVSMapper(self.bdi, cache_transcripts=True)
        self.hp = hgvs.parser.Parser()
        self.failed = []

    def test_DNAH11_hgmd(self):
        for rec in gcp_file_reader('tests/data/DNAH11-HGMD.tsv'):
            self._test_gcp_mapping(rec)

    def test_DNAH11_dbSNP(self):
        for rec in gcp_file_reader('tests/data/DNAH11-dbSNP.tsv'):
            self._test_gcp_mapping(rec)

    # TODO: Enable NEFL tests, requires adding ac=NM_006158.4
    def notest_NEFL_dbSNP(self):
        for rec in gcp_file_reader('tests/data/NEFL-dbSNP.tsv'):
            self._test_gcp_mapping(rec)

    def _test_gcp_mapping(self, rec):
        var_g = self.hp.parse_hgvs_variant(rec['HGVSg'])
        var_c = self.hp.parse_hgvs_variant(rec['HGVSc'])
        var_p = self.hp.parse_hgvs_variant(rec['HGVSp']) if rec['HGVSp'] is not None and rec['HGVSp'] != '' else None

        # g -> c
        var_c_test = self.hm.hgvsg_to_hgvsc(var_g, var_c.ac)
        self.assertEquals(str(var_c_test), str(var_c))

        # c -> g
        var_g_test = self.hm.hgvsc_to_hgvsg(var_c)
        self.assertEquals(str(var_g_test), str(var_g))

        if var_p is not None:
            # c -> p
            var_p_test = self.hm.hgvsc_to_hgvsp(var_c, var_p.ac)
            hgvs_p_test = str(var_p_test)
            self.assertEquals(hgvs_p_test, str(var_p))


if __name__ == '__main__':
    unittest.main()
