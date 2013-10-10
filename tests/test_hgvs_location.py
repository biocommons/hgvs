import unittest

import hgvs.location

class Test_Location(unittest.TestCase):
    def test_Position(self):
        pos = hgvs.location.Position(5)
        self.assertEqual( pos, 5 )
        with self.assertRaises(AssertionError):
            self.assertEqual( hgvs.location.Position(-1) )

    def test_CDSPosition(self):
        cdsp = hgvs.location.CDSPosition(5)
        self.assertEqual( cdsp.base, 5 )
        self.assertEqual( cdsp.offset, 0 )
        self.assertEqual( str(cdsp), '5' )

        cdsp = hgvs.location.CDSPosition(5,6)
        self.assertEqual( cdsp.base, 5 )
        self.assertEqual( cdsp.offset, 6 )
        self.assertEqual( str(cdsp), '5+6' )

        cdsp = hgvs.location.CDSPosition(5,-7)
        self.assertEqual( cdsp.base, 5 )
        self.assertEqual( cdsp.offset, -7 )
        self.assertEqual( str(cdsp), '5-7' )

        cdsp = hgvs.location.CDSPosition(-5,7)
        self.assertEqual( cdsp.base, -5 )
        self.assertEqual( cdsp.offset, 7 )
        self.assertEqual( str(cdsp), '-5+7' )

        cdsp = hgvs.location.CDSPosition(0,-7)
        self.assertEqual( cdsp.base, 0 )
        self.assertEqual( cdsp.offset, -7 )
        self.assertEqual( str(cdsp), '-7' )

        cdsp = hgvs.location.CDSPosition(0,+7)
        self.assertEqual( cdsp.base, 0 )
        self.assertEqual( cdsp.offset, +7 )
        self.assertEqual( str(cdsp), '*7' )


    def test_CDSInterval(self):
        cdsr = hgvs.location.Interval( hgvs.location.CDSPosition(base=12,offset=+34),
                                       hgvs.location.CDSPosition(base=56,offset=-78) )
        self.assertEqual( cdsr.start.base, 12 )
        self.assertEqual( cdsr.start.offset, 34 )
        self.assertEqual( cdsr.end.base, 56 )
        self.assertEqual( cdsr.end.offset, -78 )
        self.assertEqual( str(cdsr), '12+34_56-78' )

if __name__ == '__main__':
    unittest.main()
