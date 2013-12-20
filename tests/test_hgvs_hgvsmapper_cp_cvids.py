#
# c to p tests - all cvids
#
import os
import unittest

import test_hgvs_hgvsmapper_cp_base

class TestHgvsCToPAllCvids(test_hgvs_hgvsmapper_cp_base.TestHgvsCToPBase):

    def notest__hgvsc_to_hgvsp_cvid_plus_sample(self):
        infilename = 'cvids_all.tsv'
        outfilename = 'cvids_all.out'
        infile = os.path.join(os.path.dirname(__file__), 'data', infilename)
        outfile = os.path.join(os.path.dirname(__file__), 'data', outfilename)
        self._run_cp_test(infile, outfile)


if __name__ == '__main__':
    unittest.main()