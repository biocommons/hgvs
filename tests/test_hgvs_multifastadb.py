import os
import unittest

import hgvs.multifastadb

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
    def setUp(self):
        self.mfdb = hgvs.multifastadb.MultiFastaDB([ os.path.join(self.test_dir,s)
                                                     for s in self.sources ])


class MFDB_Test_d1f1(MFDBTBase):
    sources = ['d1/f1.fasta']

    def test_basic_file(self):
        self.assertEqual( self.mfdb.references, ['s1','s2'] )
        self.assertEqual( self.mfdb.lengths, [6,6] )

    def test_basic_sequence(self):
        self.assertEqual( self.mfdb.fetch('s1'), 'd1f1s1' )
        self.assertEqual( str(self.mfdb['s1']), 'd1f1s1' )
        self.assertEqual( self.mfdb['s1'][:2], 'd1' )


class MFDB_Test_d1(MFDBTBase):
    sources = ['d1']

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



if __name__ == '__main__':
    unittest.main()
