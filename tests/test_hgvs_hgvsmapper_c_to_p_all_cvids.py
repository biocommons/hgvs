#
# Tests for conversion of hgvs tags - real data
#
from __future__ import with_statement
import csv
import os
import re
import unittest

import hgvs.hgvsmapper as hgvsmapper
import hgvs.parser
import uta.db.transcriptdb

#INFILE = 'all_dupN.tsv'
INFILE = 'c_to_p_all_cvids_clean.tsv'
#INFILE = 'all_cvid_c_to_p_no_NONE_Met1Q_dupN.tsv'
MAPFILE = 'cvid_nm_to_np.tsv'
OUTFILE = 'all_cvid_c_to_p.out'
#OUTFILE = 'all_cvid_c_to_p_no_NONE_Met1Q_dupN_local.out'
#OUTFILE = 'all_dupN_retry.out'

class TestHgvsCToPReal(unittest.TestCase):

    _datasource = uta.db.transcriptdb.TranscriptDB()
    _mapper = hgvsmapper.HGVSMapper(_datasource, cache_transcripts=True)
    _parser = hgvs.parser.Parser()

    # def test_dbg(self):
    #     """For purposes of tesing a single result"""
    #     hgvsc = "NM_000271.4:c.3322_3323insG"
    #     hgvsp_expected = "NP_000262.2:p.Ala1108Glyfs*13"
    #     self._run_conversion(hgvsc, hgvsp_expected)

    @classmethod
    def setUpClass(cls):
        fn = os.path.join(os.path.dirname(__file__), 'data', MAPFILE)
        with open(fn, 'r') as f:
            cls._nm_np_map = {x.strip('\n').split('\t')[0]: x.strip().split('\t')[1] for x in f.readlines()}



    def notest_all_cvids(self):
        """Run all of CVID data"""
        fn = os.path.join(os.path.dirname(__file__), 'data', INFILE)
        fo = os.path.join(os.path.dirname(__file__), 'data', OUTFILE)
        ff = open(fo, 'w')
        ff.write("NM coordinates\tProtein coordinates\tConverter coordinates\tError\n")
        with open(fn, 'r') as f:
            csvreader = csv.DictReader(f, delimiter='\t')
            testdata = [(row['NM coordinates'], row['Protein coordinates']) for row in csvreader]

        regex = re.compile(r'dup[0-9]+$')
        failed_tests = []
        for x in testdata:
            print x
            (hgvsc, hgvsp_expected) = x
            hgvsc = regex.sub('dup', hgvsc) # cleanup dupN
            if not hgvsc.startswith("#"):
                try:
                    hgvsp_actual = self._run_conversion_batch(hgvsc)
                    msg = "hgvsp expected: {} actual: {}".format(hgvsp_expected, hgvsp_actual)
                    if not self._is_equivalent_hgvsp(hgvsp_expected, hgvsp_actual):
                        out = (hgvsc, hgvsp_expected, hgvsp_actual)
                        failed_tests.append(out)
                        ff.write("{}\t{}\t{}\t{}\n".format(hgvsc, hgvsp_expected, hgvsp_actual, "MISMATCH"))
                except Exception as e:
                    hgvsp_actual = "NO_OUTPUT_EXCEPTION"
                    out = (hgvsc, hgvsp_expected, hgvsp_actual)
                    failed_tests.append(out)
                    ff.write("{}\t{}\t{}\t{}\n".format(hgvsc, hgvsp_expected, hgvsp_actual, e.message))


        ff.close()
        print len(failed_tests)
        print failed_tests
        self.assertTrue(len(failed_tests) == 0)

    #
    # internal methods
    #

    def _is_equivalent_hgvsp(self, hgvsp_expected, hgvsp_actual):
        """Compare strings and account for equivalence"""
        result = True
        if hgvsp_expected != hgvsp_actual:
            if hgvsp_actual.find("Ter") != -1:
                hgvsp_actual_st = hgvsp_actual.replace("Ter", "*")
                if hgvsp_expected != hgvsp_actual_st:
                    result = False
            else:
                result = False
        return result

    def _run_conversion(self, hgvsc, hgvsp_expected):
        """Helper method to actually run the test
        :param hgvsc tag
        """
        var_c = TestHgvsCToPReal._parser.parse_hgvs_variant(hgvsc)
        var_p = TestHgvsCToPReal._mapper.hgvsc_to_hgvsp(var_c)
        #var_p.ac = self._nm_np_map[var_c.ac]
        hgvsp_actual = str(var_p)
        msg = "hgvsp expected: {} actual: {}".format(hgvsp_expected, hgvsp_actual)
        self.assertEqual(hgvsp_expected, hgvsp_actual, msg)

    def _run_conversion_batch(self, hgvsc):
        """Helper method to actually run the test
        :param hgvsc tag
        """
        var_c = TestHgvsCToPReal._parser.parse_hgvs_variant(hgvsc)
        var_p = TestHgvsCToPReal._mapper.hgvsc_to_hgvsp(var_c)
        var_p.ac = self._nm_np_map[var_c.ac]
        return str(var_p)





if __name__ == '__main__':
    unittest.main()