#
# Tests for conversion of hgvs tags - real data
#
from __future__ import with_statement
import csv
import os
import unittest

import hgvs.hgvsmapper as hgvsmapper
import hgvs.parser
import uta.db.transcriptdb

INFILE = 'eh_tests.tsv'
OUTFILE = 'eh_c_to_p.out'

class TestHgvsCToPReal(unittest.TestCase):

    _datasource = uta.db.transcriptdb.TranscriptDB()
    _mapper = hgvsmapper.HGVSMapper(_datasource, cache_transcripts=True)
    _parser = hgvs.parser.Parser()

    # def test_dbg(self):
    #     """For purposes of tesing a single result"""
    #     hgvsc = 'NM_000257.2:c.3300_3311delCGGCAGCCAGCT'
    #     hgvsp_expected = 'NP_000248.2:p.Gly1101_Leu1104del'
    #
    #     var_c = TestHgvsCToPReal._parser.parse_hgvs_variant(hgvsc)
    #     var_p = TestHgvsCToPReal._mapper.hgvsc_to_hgvsp(var_c)
    #     print "acc {}".format(var_p.ac)
    #     hgvsp_actual = str(var_p)
    #     msg = "hgvsp expected: {} actual: {}".format(hgvsp_expected, hgvsp_actual)
    #     self.assertEqual(hgvsp_expected, hgvsp_actual, msg)


    def test_eh(self):
        """Run all of Emilys data"""
        fn = os.path.join(os.path.dirname(__file__), 'data', INFILE)
        fo = os.path.join(os.path.dirname(__file__), 'data', OUTFILE)
        ff = open(fo, 'w')
        ff.write("id\tHGVSg\tHGVSc\tHGVSp\tConverterResult\tError\n")
        with open(fn, 'r') as f:
            csvreader = csv.DictReader(f, delimiter='\t')
            testdata = ((row['id'], row['HGVSg'], row['HGVSc'], row['HGVSp']) for row in csvreader)

            failed_tests = []
            for x in testdata:
                (id, hgvsg, hgvsc, hgvsp_expected) = x
                if not id.startswith("#"):
                    try:
                        var_c = TestHgvsCToPReal._parser.parse_hgvs_variant(hgvsc)
                        var_p = TestHgvsCToPReal._mapper.hgvsc_to_hgvsp(var_c)
                        hgvsp_actual = str(var_p)
                        msg = "hgvsp expected: {} actual: {}".format(hgvsp_expected, hgvsp_actual)
                        if hgvsp_expected != hgvsp_actual:
                            if hgvsp_expected.endswith("*") \
                                and hgvsp_actual.endswith("Ter") and hgvsp_expected[:-1] == hgvsp_actual[:-3]:
                                pass    # skip * vs Ter for now
                            else:
                                out = (hgvsc, hgvsp_expected, hgvsp_actual)
                                failed_tests.append(out)
                                message = "MISMATCH"
                                ff.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(id, hgvsg, hgvsc, hgvsp_expected,
                                                                           hgvsp_actual, message))
                    except Exception as e:
                        hgvsp_actual = "NO_OUTPUT_EXCEPTION"
                        out = (hgvsc, hgvsp_expected, hgvsp_actual)
                        failed_tests.append(out)
                        ff.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(id, hgvsg, hgvsc, hgvsp_expected,
                                                                   hgvsp_actual, e.message))

        ff.close()
        print len(failed_tests)
        print failed_tests
        self.assertTrue(len(failed_tests) == 0)






if __name__ == '__main__':
    unittest.main()