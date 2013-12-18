import unittest

import hgvs.location

class Test_SimplePosition(unittest.TestCase):
    def test_success(self):
        self.assertEqual( str(hgvs.location.SimplePosition(5)), '5' )
        self.assertEqual( str(hgvs.location.SimplePosition(5,uncertain=True)), '(5)' )
        self.assertEqual( str(hgvs.location.SimplePosition(None)), '?' )

    def test_failure(self):
        with self.assertRaises(AssertionError):
            self.assertEqual( hgvs.location.SimplePosition(-1), 'SHOULD FAIL' )


class Test_BaseOffsetPosition(unittest.TestCase):
    def test_success(self):
        # r.5
        cdsp = hgvs.location.BaseOffsetPosition(5)
        self.assertEqual( cdsp.datum, hgvs.location.SEQ_START )
        self.assertEqual( cdsp.base, 5 )
        self.assertEqual( cdsp.offset, 0 )
        self.assertEqual( str(cdsp), '5' )

        #r.5+6
        cdsp.offset = 6
        self.assertEqual( str(cdsp), '5+6' )

        #r.5+?
        cdsp.offset = None
        self.assertEqual( str(cdsp), '5+?' )

        #r.(5+?)
        cdsp.uncertain = True
        self.assertEqual( str(cdsp), '(5+?)' )
        
        # c.*5
        cdsp = hgvs.location.BaseOffsetPosition(5,datum=hgvs.location.CDS_END)
        self.assertEqual( cdsp.datum, hgvs.location.CDS_END )
        self.assertEqual( cdsp.base, 5 )
        self.assertEqual( cdsp.offset, 0 )
        self.assertEqual( str(cdsp), '*5' )

        cdsp.uncertain = True
        self.assertEqual( str(cdsp), '(*5)' )
        
        cdsp.offset = 7
        self.assertEqual( str(cdsp), '(*5+7)' )


class Test_AAPosition(unittest.TestCase):
    def test_AAPosition(self):
        ap = hgvs.location.AAPosition( 15, 'S' )
        self.assertEqual( ap.pos, 15 )
        self.assertEqual( str(ap), 'Ser15' )

        
class Test_Interval(unittest.TestCase):        
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
