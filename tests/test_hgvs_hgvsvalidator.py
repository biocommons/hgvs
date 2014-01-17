import unittest

import bdi.sources.uta0
import hgvs.parser
import hgvs.hgvsmapper
import hgvs.hgvsvalidator

class Test_HGVSIntrinsicValidator(unittest.TestCase):
    """Tests for internal validation"""

    def setUp(self):
        self.hp = hgvs.parser.Parser()
        self.i_valid = hgvs.hgvsvalidator.IntrinsicValidation()

    def test_check_parsing(self):
        """Test if g., r., c., p. HGVS names parse correctly"""
        self.assertTrue(self.i_valid.valid_parse('NC_000007.13:g.36561662C>T'))
        self.assertTrue(self.i_valid.valid_parse('NM_001637.3:r.1582G>A'))
        self.assertTrue(self.i_valid.valid_parse('NM_001637.3:c.1582G>A'))
        self.assertTrue(self.i_valid.valid_parse('NP_001628.1:p.Gly528Arg'))

        self.assertFalse(self.i_valid.valid_parse('NC_000007.13:z.36561662C>T'))
        self.assertFalse(self.i_valid.valid_parse('NC_000007.13:r.36561662C>T'))
        self.assertFalse(self.i_valid.valid_parse('NC_000007.13:g.36561662+2C>T'))

    def test_end_gt_start(self):
        """Test if end position is greater than start position"""
        self.assertTrue(self.i_valid.end_gt_start('NC_000007.13:g.36561662C>T'))
        self.assertTrue(self.i_valid.end_gt_start('NC_000007.13:g.36561662_36561663insT'))
        self.assertTrue(self.i_valid.end_gt_start('AC_01234.5:c.76_77insT'))
        self.assertTrue(self.i_valid.end_gt_start('AC_01234.5:c.123+54_123+55insT'))

        self.assertFalse(self.i_valid.end_gt_start('AC_01234.5:c.123+56_123+55insT'))
        self.assertFalse(self.i_valid.end_gt_start('NC_000007.13:g.36561662_36561660C>T'))

    def test_check_ins_length(self):
        """Test if ins length is 1"""
        self.assertTrue(self.i_valid.valid_ins('NC_000007.13:g.36561662_36561663insT'))
        self.assertTrue(self.i_valid.valid_ins('AC_01234.5:c.76_77insT'))
        self.assertTrue(self.i_valid.valid_ins('AC_01234.5:c.123+54_123+55insT'))
        self.assertTrue(self.i_valid.valid_ins('AC_01234.5:c.123-54_123-53insT'))

        self.assertFalse(self.i_valid.valid_ins('AC_01234.5:c.76_77delinsTT'))
        self.assertFalse(self.i_valid.valid_ins('AC_01234.5:c.123+54_123+56insT'))

    def test_check_del_length(self):
        """Test if del length agrees with position range"""
        self.assertTrue(self.i_valid.valid_del('AC_01234.5:c.76_78delACT'))
        self.assertTrue(self.i_valid.valid_del('AC_01234.5:c.123+54_123+55delTA'))  # <-- haha "delta"

        self.assertFalse(self.i_valid.valid_del('AC_01234.5:c.76_78del'))

class Test_HGVSExtrinsicValidator(unittest.TestCase):
    """Tests for external validation"""

    def setUp(self):
        self.bdi = bdi.sources.uta0.connect()
        self.hm = hgvs.hgvsmapper.HGVSMapper(self.bdi, cache_transcripts=True)
        self.hp = hgvs.parser.Parser()
        self.e_valid = hgvs.hgvsvalidator.ExtrinsicValidation()

    def test_valid_ac(self):
        """Test if accession is present in REST transcript sequence database"""
        self.assertTrue(self.e_valid.valid_ac('NM_001637.3:r.1582G>A'))

        self.assertTrue(self.e_valid.valid_ac('NM_01234.5:r.1582G>A'))



if __name__ == '__main__':
    unittest.main()

## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/invitae/hgvs)
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
