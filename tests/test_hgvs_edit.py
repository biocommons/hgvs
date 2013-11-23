import unittest

import hgvs.edit
from hgvs.exceptions import HGVSError

class Test_Edit(unittest.TestCase):
    def test_RefAlt_exceptions(self):
        with self.assertRaises(HGVSError):
            edit = str(hgvs.edit.RefAlt(None,None))

    def test_RefAlt(self):
        self.assertEqual( str(hgvs.edit.RefAlt('A','A'))	, '=' 			)
        self.assertEqual( str(hgvs.edit.RefAlt('A','T'))	, 'A>T' 		)
        self.assertEqual( str(hgvs.edit.RefAlt('AA',None))	, 'delAA' 		)
        self.assertEqual( str(hgvs.edit.RefAlt(None,'TT'))	, 'insTT' 		)
        self.assertEqual( str(hgvs.edit.RefAlt('AA','T'))	, 'delAAinsT' 	)
        self.assertEqual( str(hgvs.edit.RefAlt('A','TT'))	, 'delAinsTT' 	)

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
