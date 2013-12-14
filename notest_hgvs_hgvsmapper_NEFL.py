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

    def test_NEFL_dbSNP(self):
        tests_fn = 'tests/data/NEFL-dbSNP.tsv'
        tests_in = csv.DictReader(open(tests_fn, 'r'), delimiter='\t')
        for rec in tests_in:
            if rec['id'].startswith('#'):
                continue
            self._test_gcp_mapping(rec)
        self.assertEquals(len(self.failed), 0, self.failed)  # without this the test passes

    def _test_gcp_mapping(self, rec):
        try:
            var_g = self.hp.parse_hgvs_variant(rec['HGVSg'])
            var_c = self.hp.parse_hgvs_variant(rec['HGVSc'])
            var_p = None
            if rec['HGVSp'] != '':
                var_p = self.hp.parse_hgvs_variant(rec['HGVSp'])
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
        except Exception as e:
            error = '{}: {}'.format(rec, e.message)
            self.failed.append(error)

if __name__ == '__main__':
    unittest.main()
