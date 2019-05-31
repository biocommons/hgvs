# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import unittest

import pytest

import hgvs.edit
import hgvs.location
import hgvs.posedit


@pytest.mark.quick
@pytest.mark.models
class Test_PosEdit(unittest.TestCase):
    def test_PosEdit(self):
        pos = hgvs.location.Interval(
            hgvs.location.BaseOffsetPosition(base=12, offset=+34),
            hgvs.location.BaseOffsetPosition(base=56, offset=-78))
        edit = hgvs.edit.NARefAlt("AA", None)
        pe = hgvs.posedit.PosEdit(pos=pos, edit=edit)
        self.assertEqual(pe.format(conf={'max_ref_length': None}), "12+34_56-78delAA")


if __name__ == "__main__":
    unittest.main()

# <LICENSE>
# Copyright 2018 HGVS Contributors (https://github.com/biocommons/hgvs)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
