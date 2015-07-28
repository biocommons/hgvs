# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import unittest

from nose.plugins.attrib import attr

import hgvs.edit
import hgvs.location
from hgvs.exceptions import HGVSError


@attr(tags=["quick", "models"])
class Test_Edit(unittest.TestCase):
    def test_NARefAlt_exceptions(self):
        with self.assertRaises(HGVSError):
            edit = str(hgvs.edit.NARefAlt(None, None))

    def test_NARefAlt(self):
        self.assertEqual(str(hgvs.edit.NARefAlt('A', 'A')), '=')
        self.assertEqual(str(hgvs.edit.NARefAlt('A', 'T')), 'A>T')
        self.assertEqual(str(hgvs.edit.NARefAlt('AA', None)), 'delAA')
        self.assertEqual(str(hgvs.edit.NARefAlt(None, 'TT')), 'insTT')
        self.assertEqual(str(hgvs.edit.NARefAlt('AA', 'T')), 'delAAinsT')
        self.assertEqual(str(hgvs.edit.NARefAlt('A', 'TT')), 'delAinsTT')
        # edit types
        self.assertEqual(str(hgvs.edit.NARefAlt('A', 'A').type), 'identity')
        self.assertEqual(str(hgvs.edit.NARefAlt('A', 'T').type), 'sub')
        self.assertEqual(str(hgvs.edit.NARefAlt('AA', None).type), 'del')
        self.assertEqual(str(hgvs.edit.NARefAlt(None, 'TT').type), 'ins')
        self.assertEqual(str(hgvs.edit.NARefAlt('AA', 'T').type), 'delins')
        self.assertEqual(str(hgvs.edit.NARefAlt('A', 'TT').type), 'delins')

    def test_AARefAlt(self):
        self.assertEqual(str(hgvs.edit.AARefAlt('A', 'A')), '=')
        self.assertEqual(str(hgvs.edit.AARefAlt('A', 'T')), 'Thr')
        self.assertEqual(str(hgvs.edit.AARefAlt('AA', None)), 'del')
        self.assertEqual(str(hgvs.edit.AARefAlt(None, 'TT')), 'insThrThr')
        self.assertEqual(str(hgvs.edit.AARefAlt('', 'T')), 'delinsThr')
        self.assertEqual(str(hgvs.edit.AARefAlt('AA', 'T')), 'delinsThr')
        self.assertEqual(str(hgvs.edit.AARefAlt('A', 'TT')), 'delinsThrThr')
        # edit types
        self.assertEqual(str(hgvs.edit.AARefAlt('A', 'A').type), 'identity')
        self.assertEqual(str(hgvs.edit.AARefAlt('A', 'T').type), 'sub')
        self.assertEqual(str(hgvs.edit.AARefAlt('AA', None).type), 'del')
        self.assertEqual(str(hgvs.edit.AARefAlt(None, 'TT').type), 'ins')
        self.assertEqual(str(hgvs.edit.AARefAlt('', 'T').type), 'delins')
        self.assertEqual(str(hgvs.edit.AARefAlt('AA', 'T').type), 'delins')
        self.assertEqual(str(hgvs.edit.AARefAlt('A', 'TT').type), 'delins')

    def test_AASub(self):
        self.assertEqual(str(hgvs.edit.AASub('A', 'T')), 'Thr')
        # edit types
        self.assertEqual(str(hgvs.edit.AASub('A', 'T').type), 'sub')

    def test_AAFs(self):
        self.assertEqual(str(hgvs.edit.AAFs('A', 'T', '6')), 'ThrfsTer6')
        self.assertEqual(str(hgvs.edit.AAFs('A', 'T', '?')), 'ThrfsTer?')
        self.assertEqual(str(hgvs.edit.AAFs('A', 'T', None)), 'ThrfsTer')
        # edit types
        self.assertEqual(str(hgvs.edit.AAFs('A', 'T', '6').type), 'fs')
        self.assertEqual(str(hgvs.edit.AAFs('A', 'T', '?').type), 'fs')
        self.assertEqual(str(hgvs.edit.AAFs('A', 'T', None).type), 'fs')

    def test_AAExt(self):
        self.assertEqual(str(hgvs.edit.AAExt('A', 'V', '*', 10)), 'ValextTer10')
        self.assertEqual(str(hgvs.edit.AAExt('A', 'V', None, -10)), 'Valext-10')
        self.assertEqual(str(hgvs.edit.AAExt('A', None, None, -5)), 'ext-5')
        # edit types
        self.assertEqual(str(hgvs.edit.AAExt('A', 'V', '*', 10).type), 'ext')
        self.assertEqual(str(hgvs.edit.AAExt('A', 'V', None, -10).type), 'ext')
        self.assertEqual(str(hgvs.edit.AAExt('A', None, None, -5).type), 'ext')

    def test_Dup(self):
        self.assertEqual(str(hgvs.edit.Dup()), 'dup')
        self.assertEqual(str(hgvs.edit.Dup('T')), 'dupT')
        # edit types
        self.assertEqual(str(hgvs.edit.Dup().type), 'dup')
        self.assertEqual(str(hgvs.edit.Dup('T').type), 'dup')

    def test_Repeat(self):
        self.assertEqual(str(hgvs.edit.Repeat('CAG', 12, 34)), 'CAG(12_34)')
        # edit types
        self.assertEqual(str(hgvs.edit.Repeat('CAG', 12, 34).type), 'repeat')

    def test_Repeat_exceptions(self):
        with self.assertRaises(HGVSError):
            edit = str(hgvs.edit.Repeat('CAG', 34, 12))

    def test_Inv(self):
        self.assertEqual(str(hgvs.edit.Inv()), 'inv')
        self.assertEqual(str(hgvs.edit.Inv('AGCT')), 'invAGCT')
        self.assertEqual(str(hgvs.edit.Inv('100')), 'inv100')
        # edit types
        self.assertEqual(str(hgvs.edit.Inv().type), 'inv')

    def test_Conv(self):
        start = hgvs.location.BaseOffsetPosition(base=61, offset=-6, datum=hgvs.location.CDS_START)
        end = hgvs.location.BaseOffsetPosition(base=22, datum=hgvs.location.CDS_END)
        pos = hgvs.location.Interval(start=start, end=end)
        self.assertEqual(str(hgvs.edit.Conv('NM_001166478.1', 'c', pos)), 'conNM_001166478.1:c.61-6_*22')
        # edit types
        self.assertEqual(str(hgvs.edit.Conv('NM_001166478.1', 'c', pos).type), 'con')



if __name__ == '__main__':
    unittest.main()

## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
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
