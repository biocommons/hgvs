#
# Base class for cp tests
#
from __future__ import with_statement
import csv
import re
import unittest

import hgvs.hgvsmapper
import hgvs.parser
import uta.db.transcriptdb

def gcp_file_reader(fn):
    rdr = csv.DictReader(open(fn, 'r'), delimiter='\t')
    for rec in rdr:
        if rec['id'].startswith('#'):
            continue
        yield rec


class TestHgvsCToPBase(unittest.TestCase):

    def setUp(self):
        uta_conn = uta.db.transcriptdb.TranscriptDB()
        self._hm = hgvs.hgvsmapper.HGVSMapper(uta_conn, cache_transcripts=True)
        self._hp = hgvs.parser.Parser()
        self._failed = []

    # def test_dbg(self):
    #     """For purposes of tesing a single result"""
    #     hgvsc = 'NM_000257.2:c.3300_3311delCGGCAGCCAGCT'
    #     hgvsp_expected = 'NP_000248.2:p.Gly1101_Leu1104del'
    #     var_c = self._hp.parse_hgvs_variant(hgvsc)
    #     var_p = self._hm.hgvsc_to_hgvsp(var_c)
    #     hgvsp_actual = str(var_p)
    #     msg = "hgvsp expected: {} actual: {}".format(hgvsp_expected, hgvsp_actual)
    #     self.assertEqual(hgvsp_expected, hgvsp_actual, msg)

    #
    # internal methods
    #

    def _run_cp_test(self, infile, outfile):

        with open(outfile, 'w') as out:
            out.write("id\tHGVSg\tHGVSc\tHGVSp\tConverterResult\tError\n")
            self._dup_regex = re.compile(r'dup[0-9]+$')
            for rec in gcp_file_reader(infile):
                self._run_comparison(rec, out)
        msg = "# failed: {}".format(len(self._failed))
        self.assertTrue(len(self._failed) == 0, msg)

    def _run_comparison(self, rec, out):
        (row_id, hgvsg, hgvsc, hgvsp_expected) = (rec['id'], rec['HGVSg'], rec['HGVSc'], rec['HGVSp'])

        hgvsc = self._dup_regex.sub('dup', hgvsc) # cleanup dupN

        if not row_id.startswith("#") and hgvsc and hgvsp_expected:
            try:
                var_c = self._hp.parse_hgvs_variant(hgvsc)
                var_p = self._hm.hgvsc_to_hgvsp(var_c)
                hgvsp_actual = str(var_p)

                if hgvsp_expected != hgvsp_actual:
                    if hgvsp_expected.endswith("*") \
                        and hgvsp_actual.endswith("Ter") and hgvsp_expected[:-1] == hgvsp_actual[:-3]:
                        pass    # skip * vs Ter
                    else:
                        self._append_fail(out, row_id, hgvsg, hgvsc, hgvsp_expected, hgvsp_actual, "MISMATCH")
            except Exception as e:
                self._append_fail(out, row_id, hgvsg, hgvsc, hgvsp_expected, "NO_OUTPUT_EXCEPTION", e.message)


    def _append_fail(self, out, row_id, hgvsg, hgvsc, hgvsp_expected, hgvsp_actual, msg):
        self._failed.append((row_id, hgvsg, hgvsc, hgvsp_expected, hgvsp_actual, msg))
        out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(row_id, hgvsg, hgvsc, hgvsp_expected, hgvsp_actual, msg))

if __name__ == '__main__':
    unittest.main()