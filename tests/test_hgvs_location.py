import unittest

import hgvs.location

class Test_Location(unittest.TestCase):
    def test_Position(self):
        pos = hgvs.location.SimplePosition(5)
        self.assertEqual( pos, 5 )
        with self.assertRaises(AssertionError):
            self.assertEqual( hgvs.location.SimplePosition(-1) )

    def test_BaseOffsetPosition(self):
        # r.5
        cdsp = hgvs.location.BaseOffsetPosition(5)
        self.assertEqual( cdsp.datum, hgvs.location.SEQ_START )
        self.assertEqual( cdsp.base, 5 )
        self.assertEqual( cdsp.offset, 0 )
        self.assertEqual( str(cdsp), '5' )

        # r.5+6
        cdsp = hgvs.location.BaseOffsetPosition(5,6)
        self.assertEqual( cdsp.datum, hgvs.location.SEQ_START )
        self.assertEqual( cdsp.base, 5 )
        self.assertEqual( cdsp.offset, 6 )
        self.assertEqual( str(cdsp), '5+6' )

        # r.5-7
        cdsp = hgvs.location.BaseOffsetPosition(5,-7)
        self.assertEqual( cdsp.datum, hgvs.location.SEQ_START )
        self.assertEqual( cdsp.base, 5 )
        self.assertEqual( cdsp.offset, -7 )
        self.assertEqual( str(cdsp), '5-7' )

        # c.-5+7
        cdsp = hgvs.location.BaseOffsetPosition(-5,7,datum=hgvs.location.CDS_START)
        self.assertEqual( cdsp.datum, hgvs.location.CDS_START )
        self.assertEqual( cdsp.base, -5 )
        self.assertEqual( cdsp.offset, 7 )
        self.assertEqual( str(cdsp), '-5+7' )

        # c.*5
        cdsp = hgvs.location.BaseOffsetPosition(5,datum=hgvs.location.CDS_END)
        self.assertEqual( cdsp.datum, hgvs.location.CDS_END )
        self.assertEqual( cdsp.base, 5 )
        self.assertEqual( cdsp.offset, 0 )
        self.assertEqual( str(cdsp), '*5' )

        # c.*5+7
        cdsp = hgvs.location.BaseOffsetPosition(5,7,datum=hgvs.location.CDS_END)
        self.assertEqual( cdsp.datum, hgvs.location.CDS_END )
        self.assertEqual( cdsp.base, 5 )
        self.assertEqual( cdsp.offset, 7 )
        self.assertEqual( str(cdsp), '*5+7' )


    def Test_AAPosition(self):
        ap = hgvs.location.AAPosition( 15, 'Ser' )
        self.assertEqual( ap.pos, 15 )
        self.assertEqual( ap.aa, 'Ser' )
        self.assertEqual( str(ap), 'Ser15' )
        
    def test_Interval(self):
        ival = hgvs.location.Interval( hgvs.location.BaseOffsetPosition(base=12,offset=+34),
                                       hgvs.location.BaseOffsetPosition(base=56,offset=-78) )
        self.assertEqual( ival.start.base, 12 )
        self.assertEqual( ival.start.offset, 34 )
        self.assertEqual( ival.end.base, 56 )
        self.assertEqual( ival.end.offset, -78 )
        self.assertEqual( str(ival), '12+34_56-78' )

if __name__ == '__main__':
    unittest.main()
