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

    # def test_1(self):
    #     hgvsc = "NM_000314.4:c.68T>A"
    #     hgvsp_expected = "NP_689487.2:p.Leu23Ter"
    #     self._run_conversion(hgvsc, hgvsp_expected)
    #
    # def test_2(self):
    #     hgvsc = "NM_000314.4:c.1034T>C"
    #     hgvsp_expected = "NP_689487.2:p.Leu345Pro"
    #     self._run_conversion(hgvsc, hgvsp_expected)
    #
    # def test_3(self):
    #     hgvsc = "NM_000314.4:c.177_179delAAA"
    #     hgvsp_expected = "NP_689487.2:p.Lys60del"
    #     self._run_conversion(hgvsc, hgvsp_expected)
    #
    # def test_4(self):
    #     hgvsc = "NM_000314.4:c.966_968delAAAinsGA"
    #     hgvsp_expected = "NP_689487.2:p.Asn323Metfs*21"
    #     self._run_conversion(hgvsc, hgvsp_expected)
    #
    # def test_5(self):
    #     hgvsc = "NM_000314.4:c.824dupT"
    #     hgvsp_expected = "NP_000305.3:p.Asn276Lysfs*22"
    #     self._run_conversion(hgvsc, hgvsp_expected)
    #
    # def test_6(self):
    #     hgvsc = "NM_152274.3:c.21_22insT"
    #     hgvsp_expected = "NP_689487.2:p.Gly8Trpfs*50"
    #     self._run_conversion(hgvsc, hgvsp_expected)

    # def test_dbg(self):
    #     hgvsc = "NM_000314.4:c.987_990delTAAA"
    #     hgvsp_expected = "NP_000305.3:p.Asn329Lysfs*14"
    #     self._run_conversion(hgvsc, hgvsp_expected)


    def test_eh(self):
        """Run all of Emilys data"""
        fn = os.path.join(os.path.dirname(__file__), 'data', 'hgvs_c_to_p_eh.tsv')
        fo = os.path.join(os.path.dirname(__file__), 'data', 'hgvs_c_to_p_eh.out')
        ff = open(fo, 'w')
        with open(fn, 'r') as f:
            dialect = csv.Sniffer().sniff(f.read(1024))
            f.seek(0)
            csvreader = csv.DictReader(f, dialect=dialect)
            testdata = [(row['NM'], row['NP']) for row in csvreader]

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