import unittest

import hgvs.hgvsposition

class Test_HGVSPosition(unittest.TestCase):
    def test_hgvsposition(self):
        var = hgvs.hgvsposition.HGVSPosition(
            seqref='NM_01234.5',
            type='c',
            pos=hgvs.location.Interval( hgvs.location.BaseOffsetPosition(base=12,offset=+34),
                                        hgvs.location.BaseOffsetPosition(base=56,offset=-78) ) )

        self.assertEqual( str(var) , 'NM_01234.5:c.12+34_56-78' )

if __name__ == '__main__':
    unittest.main()
