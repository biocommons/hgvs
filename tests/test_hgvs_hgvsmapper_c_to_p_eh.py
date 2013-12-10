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

# map to handle stopgap
PROTEIN_ACCESSION_MAP = {
    "MD5_ef02ef58": "NP_689487.2",
    "MD5_bdd27960": "NP_079016.2",
    "MD5_51f440ba": "NP_006149.2",
    "MD5_6009ae51": "NP_000390.2",
    "MD5_7e4eaf25": "NP_000248.2",
    "MD5_ce4981e5": "NP_000305.3",
    "MD5_5622254c": "NP_003051.1",
}


class TestHgvsCToPReal(unittest.TestCase):

    _datasource = uta.db.transcriptdb.TranscriptDB()
    _mapper = hgvsmapper.HGVSMapper(_datasource, cache_transcripts=True)
    _parser = hgvs.parser.Parser()

    def test_dbg(self):
        """For purposes of tesing a single result"""
        hgvsc = "NM_000314.4:c.800dupA"
        hgvsp_expected = "NP_000305.3:p.Asp268Glyfs*9"
        self._run_conversion(hgvsc, hgvsp_expected)


    def test_eh(self):
        """Run all of Emilys data"""
        fn = os.path.join(os.path.dirname(__file__), 'data', 'eh_tests.tsv')
        fo = os.path.join(os.path.dirname(__file__), 'data', 'hgvs_c_to_p_eh.out')
        ff = open(fo, 'w')
        with open(fn, 'r') as f:
            csvreader = csv.DictReader(f, delimiter='\t')
            testdata = [(row['NM coordinates'], row['Protein coordinates']) for row in csvreader]

        failed_tests = []
        for x in testdata:
            (hgvsc, hgvsp_expected) = x
            try:
                hgvsp_actual = self._run_conversion_batch(hgvsc)
                msg = "hgvsp expected: {} actual: {}".format(hgvsp_expected, hgvsp_actual)
                if hgvsp_expected != hgvsp_actual:
                    if hgvsp_expected.endswith("*") \
                        and hgvsp_actual.endswith("Ter") and hgvsp_expected[:-1] == hgvsp_actual[:-3]:
                        pass
                    else:
                        out = (hgvsc, hgvsp_expected, hgvsp_actual)
                        failed_tests.append(out)
                        ff.write("{}\t{}\t{}\n".format(hgvsc, hgvsp_expected, hgvsp_actual))
            except Exception as e:
                print e.message
                print "exception processing {}".format(x)
                hgvsp_actual = "NO_OUTPUT_EXCEPTION"
                ff.write("{}\t{}\t{}\n".format(hgvsc, hgvsp_expected, hgvsp_actual))

        ff.close()
        print len(failed_tests)
        print failed_tests
        self.assertTrue(len(failed_tests) == 0)

    #
    # internal methods
    #

    def _run_conversion(self, hgvsc, hgvsp_expected):
        """Helper method to actually run the test
        :param hgvsc tag
        """
        var_c = TestHgvsCToPReal._parser.parse_hgvs_variant(hgvsc)
        var_p = TestHgvsCToPReal._mapper.hgvsc_to_hgvsp(var_c)
        var_p.ac = PROTEIN_ACCESSION_MAP[var_p.ac]
        hgvsp_actual = str(var_p)
        msg = "hgvsp expected: {} actual: {}".format(hgvsp_expected, hgvsp_actual)
        self.assertEqual(hgvsp_expected, hgvsp_actual, msg)




    def _run_conversion_batch(self, hgvsc):
        """Helper method to actually run the test
        :param hgvsc tag
        """
        var_c = TestHgvsCToPReal._parser.parse_hgvs_variant(hgvsc)
        var_p = TestHgvsCToPReal._mapper.hgvsc_to_hgvsp(var_c)
        var_p.ac = PROTEIN_ACCESSION_MAP[var_p.ac]
        return str(var_p)





if __name__ == '__main__':
    unittest.main()