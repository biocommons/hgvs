import unittest

from nose.plugins.attrib import attr

from hgvs.dataproviders.multifastadb import MultiFastaDB
from hgvs.exceptions import HGVSValidationError
import hgvs.dataproviders.uta
import hgvs.variantmapper
import hgvs.parser
import hgvs.validator

db_dir = ['tests/data/sample_data']
hdp = hgvs.dataproviders.uta.connect()
mfdb = MultiFastaDB(db_dir, use_meta_index=True)


class Test_HGVSValidator(unittest.TestCase):
    """Validator wrapper class tests (most testing is handled by the component classes)"""

    def setUp(self):
        self.hp = hgvs.parser.Parser()
        self.vr = hgvs.validator.Validator(hdp,mfdb)

    def test_wrapper(self):
        """Test that validator wrapper is working"""
        self.assertTrue(self.vr.validate(self.hp.parse_hgvs_variant('NM_001005405.2:c.6C>A')))


@attr(tags=["quick","validation"])
class Test_HGVSIntrinsicValidator(unittest.TestCase):
    """Tests for internal validation"""

    def setUp(self):
        self.hp = hgvs.parser.Parser()
        self.validate_int = hgvs.validator.IntrinsicValidator()

    def test_start_lte_end(self):
        """Test if start position is less <= end position"""
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('NC_000007.13:g.36561662C>T')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('NC_000007.13:g.36561662_36561663insT')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('NM_01234.5:c.76_77insT')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123+54_123+55insT')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123+54A>T')))

        with self.assertRaisesRegexp(HGVSValidationError, hgvs.validator.BASE_RANGE_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('NC_000007.13:g.36561664_36561663A>T'))

        with self.assertRaisesRegexp(HGVSValidationError, hgvs.validator.BASE_RANGE_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('NM_000277.1:c.*1_2delAG'))

        with self.assertRaisesRegexp(HGVSValidationError, hgvs.validator.OFFSET_RANGE_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123+56_123+55A>T'))

    def test_ins_length_is_one(self):
        """Test if insertion length = 1"""
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('NC_000007.13:g.36561662_36561663insT')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.76_77insT')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123+54_123+55insT')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123-54_123-53insT')))

        with self.assertRaisesRegexp(HGVSValidationError, hgvs.validator.INS_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.76_78insTT'))

        with self.assertRaisesRegexp(HGVSValidationError, hgvs.validator.INS_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123+54_123+56insT'))

    def test_del_length(self):
        """Test if del length agrees with position range"""
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.76_78delACT')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123+54_123+55delTA')))  # <-- haha "delta"

        with self.assertRaisesRegexp(HGVSValidationError, hgvs.validator.DEL_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.76_78del'))

        with self.assertRaisesRegexp(HGVSValidationError, hgvs.validator.DEL_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.76_78delACTACAT'))

    def test_sub(self):
        """Test substitution ref != alt"""
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('NC_000007.13:g.36561662C>T')))
        self.assertTrue(self.validate_int.validate(self.hp.parse_hgvs_variant('AC_01234.5:c.123+54A>T')))

        with self.assertRaisesRegexp(HGVSValidationError, hgvs.validator.SUB_ERROR_MSG):
            self.validate_int.validate(self.hp.parse_hgvs_variant('NC_000007.13:g.36561662_36561663T>T'))


@attr(tags=["validation"])
class Test_HGVSExtrinsicValidator(unittest.TestCase):
    """Tests for external validation"""

    def setUp(self):
        self.hp = hgvs.parser.Parser()
        self.validate_ext = hgvs.validator.ExtrinsicValidator(hdp,mfdb)

    def test_valid_ref(self):
        """Test if reference seqeuence is valid. Uses sample_data in tests directory"""
        self.assertTrue(self.validate_ext.validate(self.hp.parse_hgvs_variant('NM_001005405.2:c.6C>A')))
        self.assertTrue(self.validate_ext.validate(self.hp.parse_hgvs_variant('NM_001005405.2:c.-38T>A')))
        self.assertTrue(self.validate_ext.validate(self.hp.parse_hgvs_variant('NM_001005405.2:c.*3C>G')))
        self.assertTrue(self.validate_ext.validate(self.hp.parse_hgvs_variant('NM_001005405.2:c.435_440delCTGCTG')))
        #self.assertTrue(self.validate_ext.validate(self.hp.parse_hgvs_variant('NP_001005405.1:p.Gly2Ser')))

        with self.assertRaisesRegexp(HGVSValidationError, hgvs.validator.SEQ_ERROR_MSG):
            self.validate_ext.validate(self.hp.parse_hgvs_variant('NM_001005405.2:c.435_440delCTGCT'))


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
