import unittest

import hgvs.edit
from hgvs.exceptions import HGVSError

class Test_Edit(unittest.TestCase):
    def test_NARefAlt_exceptions(self):
        with self.assertRaises(HGVSError):
            edit = str(hgvs.edit.NARefAlt(None,None))

    def test_NARefAlt(self):
        self.assertEqual( str(hgvs.edit.NARefAlt('A','A'  ))  , '='             )
        self.assertEqual( str(hgvs.edit.NARefAlt('A','T'  ))  , 'A>T'           )
        self.assertEqual( str(hgvs.edit.NARefAlt('AA',None))  , 'delAA'         )
        self.assertEqual( str(hgvs.edit.NARefAlt(None,'TT'))  , 'insTT'         )
        self.assertEqual( str(hgvs.edit.NARefAlt('AA','T' ))  , 'delAAinsT'     )
        self.assertEqual( str(hgvs.edit.NARefAlt('A','TT' ))  , 'delAinsTT'     )


    def test_AARefAlt(self):
        self.assertEqual( str(hgvs.edit.AARefAlt('A','A'  ))  , '='             )
        self.assertEqual( str(hgvs.edit.AARefAlt('A','T'  ))  , 'Thr'           )
        self.assertEqual( str(hgvs.edit.AARefAlt('AA',None))  , 'del'           )
        self.assertEqual( str(hgvs.edit.AARefAlt(None,'TT'))  , 'insThrThr'     )
        self.assertEqual( str(hgvs.edit.AARefAlt('AA','T' ))  , 'delinsThr'     )
        self.assertEqual( str(hgvs.edit.AARefAlt('A','TT' ))  , 'delinsThrThr'  )

    def test_Dup(self):
        self.assertEqual( str(hgvs.edit.Dup())				, 'dup' 		)
        self.assertEqual( str(hgvs.edit.Dup('T'))			, 'dupT' 		)
        
    def test_Repeat(self):
        self.assertEqual( str(hgvs.edit.Repeat('CAG',12,34)), 'CAG(12_34)' 	)

    def test_Repeat_exceptions(self):
        with self.assertRaises(HGVSError):
            edit = str(hgvs.edit.Repeat('CAG',34,12))

if __name__ == '__main__':
    unittest.main()
