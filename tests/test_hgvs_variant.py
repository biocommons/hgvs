import unittest

import hgvs.variant

class Test_Variant(unittest.TestCase):
    def test_variant(self):
        var = hgvs.variant.SequenceVariant(ac='AC',type='B',posedit='1234DE>FG')
        self.assertEqual( str(var) , 'AC:B.1234DE>FG' )

if __name__ == '__main__':
    unittest.main()
