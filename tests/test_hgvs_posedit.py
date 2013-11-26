import unittest

import hgvs.edit
import hgvs.location
import hgvs.posedit

class Test_PosEdit(unittest.TestCase):
    def test_PosEdit(self):
        pos = hgvs.location.Interval(
                hgvs.location.BaseOffsetPosition(base=12,offset=+34),
                hgvs.location.BaseOffsetPosition(base=56,offset=-78) )
        edit = hgvs.edit.NARefAlt('AA',None)
        pe = hgvs.posedit.PosEdit( pos=pos, edit=edit )
        self.assertEqual( str(pe), '12+34_56-78delAA' )

if __name__ == '__main__':
    unittest.main()
