import unittest

import bdi.sources.uta0
import hgvs.parser
import hgvs.hgvsmapper
import hgvs.validator

class Test_HGVSIntrinsicValidator(unittest.TestCase):
    """Tests for internal validation"""

    def setUp(self):
        self.hp = hgvs.parser.Parser()
        self.validate_int = hgvs.validator.IntrinsicValidation()

    def test_start_lte_end(self):
        """Test if start position is less <= end position"""
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('NC_000007.13:g.36561662C>T')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('NC_000007.13:g.36561662_36561663insT')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('NM_01234.5:c.76_77insT')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123+54_123+55insT')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123+54A>T')))

        with self.assertRaisesRegexp(Exception, self.validate_int.BASE_RANGE_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('NC_000007.13:g.36561664_36561663A>T'))

        with self.assertRaisesRegexp(Exception, self.validate_int.BASE_RANGE_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('NM_000277.1:c.*1_2delAG'))

        with self.assertRaisesRegexp(Exception, self.validate_int.OFFSET_RANGE_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123+56_123+55A>T'))

    def test_ins_length_is_one(self):
        """Test if insertion length = 1"""
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('NC_000007.13:g.36561662_36561663insT')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.76_77insT')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123+54_123+55insT')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123-54_123-53insT')))

        with self.assertRaisesRegexp(Exception, self.validate_int.INS_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.76_78insTT'))

        with self.assertRaisesRegexp(Exception, self.validate_int.INS_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123+54_123+56insT'))

    def test_del_length(self):
        """Test if del length agrees with position range"""
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.76_78delACT')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123+54_123+55delTA')))  # <-- haha "delta"

        with self.assertRaisesRegexp(Exception, self.validate_int.DEL_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.76_78del'))

        with self.assertRaisesRegexp(Exception, self.validate_int.DEL_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.76_78delACTACAT'))

#class Test_HGVSExtrinsicValidator(unittest.TestCase):
#    """Tests for external validation"""
#
#    def setUp(self):
#        self.bdi = bdi.sources.uta0.connect()
#        self.hm = hgvs.hgvsmapper.HGVSMapper(self.bdi, cache_transcripts=True)
#        self.hp = hgvs.parser.Parser()
#        self.e_valid = hgvs.hgvsvalidator.ExtrinsicValidation()
#
#    def test_valid_ac(self):
#        """Test if accession is present in REST transcript sequence database"""
#        self.assertTrue(self.e_valid.valid_ac('NM_001637.3:r.1582G>A'))
#
#        self.assertTrue(self.e_valid.valid_ac('NM_01234.5:r.1582G>A'))



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
