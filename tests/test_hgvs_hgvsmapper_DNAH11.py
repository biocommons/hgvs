import unittest, csv
import hgvs.hgvsmapper as hgvsmapper
import hgvs.parser
import uta.db.transcriptdb


class Test_HGVSMapper(unittest.TestCase):
    def setUp(self):
        uta_conn = uta.db.transcriptdb.TranscriptDB()
        self.hm = hgvs.hgvsmapper.HGVSMapper(uta_conn, cache_transcripts=True)
        self.hp = hgvs.parser.Parser()
        self.failed = []

    def test_DNAH11_hgmd(self):
        tests_fn = 'tests/data/DNAH11-HGMD.tsv'
        tests_in = csv.DictReader(open(tests_fn, 'r'), delimiter='\t')
        for rec in tests_in:
            if rec['id'].startswith('#'):
                continue
            var_g = self.hp.parse_hgvs_variant(rec['HGVSg'])
            var_c = self.hp.parse_hgvs_variant(rec['HGVSc'])
            if rec['HGVSp'] == 'x':
                var_p = self.hp.parse_hgvs_variant(rec['HGVSp'])
                self._test_gcp_mapping(var_g, var_c, var_p)
            else:
                self._test_gcp_mapping(var_g, var_c)
        print self.failed

    def test_DNAH11_dbSNP(self):
        tests_fn = 'tests/data/DNAH11-dbSNP.tsv'
        tests_in = csv.DictReader(open(tests_fn, 'r'), delimiter='\t')
        for rec in tests_in:
            if rec['id'].startswith('#'):
                continue
            var_g = self.hp.parse_hgvs_variant(rec['HGVSg'])
            var_c = self.hp.parse_hgvs_variant(rec['HGVSc'])
            if rec['HGVSp'] == 'x':
                var_p = self.hp.parse_hgvs_variant(rec['HGVSp'])
                self._test_gcp_mapping(var_g, var_c, var_p)
            else:
                self._test_gcp_mapping(var_g, var_c)
        print self.failed


    def _test_gcp_mapping(self, var_g, var_c, var_p=None):
        try:
            # g -> c
            var_c_test = self.hm.hgvsg_to_hgvsc(var_g, var_c.ac)
            self.assertEquals(str(var_c_test), str(var_c))
            # c -> g
            var_g_test = self.hm.hgvsc_to_hgvsg(var_c)
            self.assertEquals(str(var_g_test), str(var_g))
            if var_p is not None:
                # g -> p
                var_p_test = self.hm.hgvsc_to_hgvsp(self.hm.hgvsg_to_hgvsc(var_g, var_c.ac))
                self.assertEquals(str(var_p_test), str(var_p))
        except Exception, e:
            self.failed.append(e)

if __name__ == '__main__':
    unittest.main()

