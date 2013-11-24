import unittest

import hgvs.utils

from hgvs.exceptions import HGVSError

class Test_Utils(unittest.TestCase):
    all_aa1 = 'ACDEFGHIKLMNPQRSTVWY'
    all_aa3 = 'AlaArgAsnAspCysGlnGluGlyHisIleLeuLysMetPheProSerThrTrpTyrVal'

    aa1_seq = 'YWVTSRQPNMLKIHGFEDCA'
    aa3_seq = 'TyrTrpValThrSerArgGlnProAsnMetLeuLysIleHisGlyPheGluAspCysAla'


    def test_aa1_lut(self):
        self.assertEqual( len(hgvs.utils.aa1_to_aa3_lut), 20 )
        self.assertEqual( ''.join(sorted(hgvs.utils.aa1_to_aa3_lut.keys())),   self.all_aa1 )
        self.assertEqual( ''.join(sorted(hgvs.utils.aa1_to_aa3_lut.values())), self.all_aa3 )

    def test_aa3_lut(self):
        self.assertEqual( len(hgvs.utils.aa3_to_aa1_lut), 20 )
        self.assertEqual( ''.join(sorted(hgvs.utils.aa3_to_aa1_lut.values())), self.all_aa1 )
        self.assertEqual( ''.join(sorted(hgvs.utils.aa3_to_aa1_lut.keys())),   self.all_aa3 )


    def test_aa1_to_aa3(self):
        self.assertEqual( hgvs.utils.aa3_to_aa1(self.aa3_seq), self.aa1_seq )

    def test_aa3_to_aa1(self):
        self.assertEqual( hgvs.utils.aa1_to_aa3(self.aa1_seq), self.aa3_seq )


    def test_aa_to_aa1(self):
        self.assertEqual( hgvs.utils.aa_to_aa1('')		, ''   )
        self.assertEqual( hgvs.utils.aa_to_aa1('A')     , 'A'  )
        self.assertEqual( hgvs.utils.aa_to_aa1('AC')    , 'AC' )
        self.assertEqual( hgvs.utils.aa_to_aa1('Ala')   , 'A'  )
        self.assertEqual( hgvs.utils.aa_to_aa1('AlaCys'), 'AC' )
        self.assertEqual( hgvs.utils.aa_to_aa1(self.aa1_seq), self.aa1_seq )
        self.assertEqual( hgvs.utils.aa_to_aa1(self.aa3_seq), self.aa1_seq )

    def test_aa_to_aa3(self):
        self.assertEqual( hgvs.utils.aa_to_aa3('')      , ''       )
        self.assertEqual( hgvs.utils.aa_to_aa3('A')     , 'Ala'    )
        self.assertEqual( hgvs.utils.aa_to_aa3('AC')    , 'AlaCys' )
        self.assertEqual( hgvs.utils.aa_to_aa3('Ala')   , 'Ala'    )
        self.assertEqual( hgvs.utils.aa_to_aa3('AlaCys'), 'AlaCys' )
        self.assertEqual( hgvs.utils.aa_to_aa3(self.aa1_seq), self.aa3_seq )
        self.assertEqual( hgvs.utils.aa_to_aa3(self.aa3_seq), self.aa3_seq )


if __name__ == '__main__':
    unittest.main()
