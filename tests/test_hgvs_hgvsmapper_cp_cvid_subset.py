#
# c to p tests - cvid subset plus some additional tests (from Emily Hare)
#
import os
import unittest

import test_hgvs_hgvsmapper_cp_base

class TestHgvsCToPSubsetCvidsPlus(test_hgvs_hgvsmapper_cp_base.TestHgvsCToPBase):

    def test__hgvsc_to_hgvsp_cvid_plus_sample(self):
        infilename = 'cvid_subset_plus.tsv'
        outfilename = 'cvid_subset_plus.out'
        infile = os.path.join(os.path.dirname(__file__), 'data', infilename)
        outfile = os.path.join(os.path.dirname(__file__), 'data', outfilename)
        self._run_cp_test(infile, outfile)


if __name__ == '__main__':
    unittest.main()