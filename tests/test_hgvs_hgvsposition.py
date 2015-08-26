# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import unittest

from nose.plugins.attrib import attr

import hgvs.hgvsposition


class Test_HGVSPosition(unittest.TestCase):
    @attr(tags=["quick", "models"])
    def test_hgvsposition(self):
        var = hgvs.hgvsposition.HGVSPosition(
            ac='NM_01234.5',
            type='c',
            pos=hgvs.location.Interval(hgvs.location.BaseOffsetPosition(base=12,
                                                                        offset=+34),
                                       hgvs.location.BaseOffsetPosition(base=56,
                                                                        offset=-78)))

        self.assertEqual(str(var), 'NM_01234.5:c.12+34_56-78')


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
