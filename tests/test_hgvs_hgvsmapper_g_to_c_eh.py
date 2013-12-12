#
# Tests for conversion of hgvs tags
#
from __future__ import with_statement
import csv
import os
import unittest

import hgvs.hgvsmapper as hgvsmapper
import hgvs.parser
import uta.db.transcriptdb

INFILE = 'eh_tests.tsv'


class TestHgvsGToCReal(unittest.TestCase):

    _datasource = uta.db.transcriptdb.TranscriptDB()
    _mapper = hgvsmapper.HGVSMapper(_datasource, cache_transcripts=True)
    _parser = hgvs.parser.Parser()


    def notest_eh_g_to_c(self):
        """Run all of Emilys data - g to c"""
        self._run_batch(TestHgvsGToCReal._mapper.hgvsg_to_hgvsc, "NC coordinates", "NM coordinates", "hgvs_g_to_c.out")

    #
    # internal methods
    #

    def _run_batch(self, func, input_header, expected_header, outfilename):
        """Run all of Emilys data"""
        fn = os.path.join(os.path.dirname(__file__), 'data', INFILE)
        fo = os.path.join(os.path.dirname(__file__), 'data', outfilename)
        ff = open(fo, 'w')
        with open(fn, 'r') as f:
            csvreader = csv.DictReader(f, delimiter='\t')
            testdata = [(row[input_header], row[expected_header]) for row in csvreader]

        failed_tests = []
        for x in testdata:
            (hgvs_in, hgvs_expected) = x
            try:
                var_g = TestHgvsGToCReal._parser.parse_hgvs_variant(hgvs_in)
                var_c_expected = TestHgvsGToCReal._parser.parse_hgvs_variant(hgvs_expected)
                hgvs_actual = str(func(var_g, var_c_expected.ac))
                msg = "{}\t{}\t{}\n".format(hgvs_in, hgvs_expected, hgvs_actual)
                if hgvs_expected != hgvs_actual:
                    out = (hgvs_in, hgvs_expected, hgvs_actual)
                    failed_tests.append(out)
                    ff.write(msg)
            except Exception as e:
                print e.message
                print "exception processing {}".format(x)
                hgvs_actual = "NO_OUTPUT_EXCEPTION"
                ff.write("{}\t{}\t{}\n".format(hgvs_in, hgvs_expected, hgvs_actual))

        ff.close()
        print len(failed_tests)
        print failed_tests
        self.assertTrue(len(failed_tests) == 0)




if __name__ == '__main__':
    unittest.main()