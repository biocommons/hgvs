import os
import unittest

import hgvs.dataproviders.multifastadb

# ==> d1/f1.fasta <==
# >s1
# d1f1s1
# >s2
# d1f1s2
# 
# ==> d1/f2.fasta <==
# >s3
# d1f2s3
# >s4
# d1f2s4
# 
# ==> d2/f3.fa <==
# >s1
# d2f3s1
# >s5
# d2f3s5
# 
# ==> d2/f4.fa <==
# >s3
# d2f4s3
# >s6
# d2f4s6


class MFDBTBase(unittest.TestCase):
    test_dir = os.path.join(os.path.dirname(__file__),'data','multifastadb')
    print test_dir
    def setUp(self):
        self.mfdb = hgvs.dataproviders.multifastadb.MultiFastaDB([ os.path.join(self.test_dir,s)
                                                     for s in self.sources ], use_meta_index=self.use_meta_index)


class MFDB_Test_d1f1(MFDBTBase):
    sources = ['d1/f1.fasta']
    use_meta_index = False

    def test_basic_file(self):
        self.assertEqual( self.mfdb.references, ['s1','s2'] )
        self.assertEqual( self.mfdb.lengths, [6,6] )

    def test_basic_sequence(self):
        self.assertEqual( self.mfdb.fetch('s1'), 'd1f1s1' )
        self.assertEqual( str(self.mfdb['s1']), 'd1f1s1' )
        self.assertEqual( self.mfdb['s1'][:2], 'd1' )


class MFDB_Test_d1(MFDBTBase):
    sources = ['d1']
    use_meta_index = False

    def test_basic_file(self):
        self.assertEqual( self.mfdb.references, ['s1','s2','s3','s4'] )
        self.assertEqual( self.mfdb.lengths, [6,6,6,6] )

    def test_basic_sequence(self):
        self.assertEqual( self.mfdb.fetch('s1'), 'd1f1s1' )
        self.assertEqual( str(self.mfdb['s1']), 'd1f1s1' )
        self.assertEqual( self.mfdb['s1'][:2], 'd1' )
        self.assertEqual( self.mfdb.fetch('s3'), 'd1f2s3' )
        self.assertEqual( str(self.mfdb['s3']), 'd1f2s3' )
        self.assertEqual( self.mfdb['s3'][:2], 'd1' )

class MFDB_Test_d2d1(MFDBTBase):
    sources = ['d2','d1']
    use_meta_index = False

    def test_basic_file(self):
        self.assertEqual( self.mfdb.references, ['s1','s5','s3','s6','s1','s2','s3','s4'] )
        self.assertEqual( self.mfdb.lengths, [6,6,6,6,6,6,6,6] )

    def test_basic_sequence(self):
        # s1 appears twice; should get d2f3 version
        self.assertEqual( self.mfdb.fetch('s1'), 'd2f3s1' )
        self.assertEqual( str(self.mfdb['s1']), 'd2f3s1' )
        self.assertEqual( self.mfdb['s1'][:2], 'd2' )

        self.assertEqual( self.mfdb.fetch('s5'), 'd2f3s5' )
        self.assertEqual( str(self.mfdb['s5']), 'd2f3s5' )
        self.assertEqual( self.mfdb['s5'][:2], 'd2' )

class MFDB_Test_Actual(MFDBTBase):
    sources = ['real_examples']
    use_meta_index = True

    def test_basic_file(self):
        self.assertEquals( self.mfdb.references, ['gi|53292629|ref|NP_001005405.1|',
                                                  'gi|52317162|ref|NP_001004713.1|',
                                                  'gi|123173798|ref|NM_001005405.2|',
                                                  'gi|52317161|ref|NM_001004713.1|',
                                                  'gi|296923737|ref|NG_021245.2|',
                                                  'gi|226510210|ref|NG_011805.1|',
                                                  ])
        self.assertEquals( self.mfdb.lengths, [156,355,1029,1068,96420,1169911] )

    def test_basic_sequence(self):
        self.assertEquals( self.mfdb.fetch('gi|53292629|ref|NP_001005405.1|'), 'MGCCGCSGGCGSGCGGCGSGSGGCGSGCGGCGSSCCVPICCCKPVCCCVPACSCSSCGSCGGSKGGCGSCGSSKGGCGSCGCSQSNCCKPCCSSSGCGSFCCQSSCSKPCCCQSSCCQSSCCKPCCCQSSCCQSSCFKPCCCQSSCCVPVCCQCKI')
        self.assertEquals( self.mfdb.fetch('NP_001005405.1'), 'MGCCGCSGGCGSGCGGCGSGSGGCGSGCGGCGSSCCVPICCCKPVCCCVPACSCSSCGSCGGSKGGCGSCGSSKGGCGSCGCSQSNCCKPCCSSSGCGSFCCQSSCSKPCCCQSSCCQSSCCKPCCCQSSCCQSSCFKPCCCQSSCCVPVCCQCKI')
        self.assertEquals( self.mfdb.fetch('NP_001005405.1',0,3), 'MGC')

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
