import unittest

import hgvs.location

class Test_Position(unittest.TestCase):
    def test_Position(self):
        self.assertEqual( hgvs.location.Position(5), 5 )
        with self.assertRaises(AssertionError):
            self.assertEqual( hgvs.location.Position(-1) )

    def test_CDSPosition(self):
        cdsp = hgvs.location.CDSPosition(5,10)
        self.assertEqual( cdsp.base, 5 )
        self.assertEqual( cdsp.offset, 10 )

    def test_Range(self):
        cdsr = hgvs.location.Range( hgvs.location.CDSPosition(base=55,offset=-10),
                                    hgvs.location.CDSPosition(base=56,offset=-8) )
        self.assertEqual( cdsr.start.base, 55 )
        self.assertEqual( cdsr.start.offset, -10 )
        self.assertEqual( cdsr.end.base, 56 )
        self.assertEqual( cdsr.end.offset, -8 )

if __name__ == '__main__':
    unittest.main()
