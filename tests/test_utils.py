import unittest

import hgvs.utils

from hgvs.exceptions import HGVSError

class Test_Utils(unittest.TestCase):
    # set of aa1 and aa3 -- must be in alpha order
    all_aa1 = '*ACDEFGHIKLMNPQRSTUVWXY'
    all_aa3 = 'AlaArgAsnAspCysGlnGluGlyHisIleLeuLysMetPheProSecSerTerThrTrpTyrValXaa'

    # sequence of aa1 and aa3 in corresponding order
    aa1_seq = 'YWVTSRQPNMLKIHGFEDCAXU*'
    aa3_seq = 'TyrTrpValThrSerArgGlnProAsnMetLeuLysIleHisGlyPheGluAspCysAlaXaaSecTer'


    def test_aa1_lut(self):
        self.assertEqual( len(hgvs.utils.aa1_to_aa3_lut), 23 )
        self.assertEqual( ''.join(sorted(hgvs.utils.aa1_to_aa3_lut.keys())),   self.all_aa1 )
        self.assertEqual( ''.join(sorted(hgvs.utils.aa1_to_aa3_lut.values())), self.all_aa3 )

    def test_aa3_lut(self):
        self.assertEqual( len(hgvs.utils.aa3_to_aa1_lut), 23 )
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

## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/invitae/hgvs)
## 
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
## 
##     http://www.apache.org/licenses/LICENSE-2.0
## 
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
