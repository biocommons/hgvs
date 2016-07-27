# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import unittest

from nose.plugins.attrib import attr

import hgvs.exceptions
import hgvs.intervalmapper


@attr(tags=["quick"])
class Test_IntervalMapper(unittest.TestCase):
    longMessage = True

    def test_interval_valid(self):
        interval = hgvs.intervalmapper.Interval(3, 4)
        st = str(interval)
        self.assertEqual(1, interval.len)
        self.assertEqual(interval.__repr__(), 'Interval(start_i=3,end_i=4)')

    def test_interval_invalid_start_gt_end(self):
        with self.assertRaises(hgvs.exceptions.HGVSInvalidIntervalError):
            interval = hgvs.intervalmapper.Interval(4, 3)

    def test_intervalpair_valid_match(self):
        self._check_valid_intervalpair(3, 4, 6, 7)

    def test_intervalpair_valid_insertion(self):
        self._check_valid_intervalpair(3, 3, 6, 7)

    def test_intervalpair_valid_deletion(self):
        self._check_valid_intervalpair(3, 4, 6, 6)

    def test_intervalpair_invalid(self):
        iv1 = hgvs.intervalmapper.Interval(3, 4)
        iv2 = hgvs.intervalmapper.Interval(6, 9)
        with self.assertRaises(hgvs.exceptions.HGVSInvalidIntervalError):
            intervalpair = hgvs.intervalmapper.IntervalPair(iv1, iv2)

    def test_intervalmapper_valid_ranges(self):
        ivm = self._build_mock_intervalmapper()

        self.assertEqual((9, 8), (ivm.ref_len, ivm.tgt_len))
        self.assertEqual((1, 8), ivm.map_tgt_to_ref(1, 7))
        self.assertEqual((1, 6), ivm.map_ref_to_tgt(1, 7))

    def test_intervalmapper_invalid_start_gt_end(self):
        ivm = self._build_mock_intervalmapper()
        with self.assertRaises(AssertionError):
            (s, e) = ivm.map_tgt_to_ref(5, 4)

    def test_intervalmapper_invalid_out_of_range(self):
        ivm = self._build_mock_intervalmapper()
        with self.assertRaises(hgvs.exceptions.HGVSInvalidIntervalError):
            (s, e) = ivm.map_tgt_to_ref(1, 200)
        with self.assertRaises(hgvs.exceptions.HGVSInvalidIntervalError):
            (s, e) = ivm.map_tgt_to_ref(0, 7)

    #
    # internal methods
    #

    def _check_valid_intervalpair(self, s1, e1, s2, e2):
        iv1 = hgvs.intervalmapper.Interval(s1, e1)
        iv2 = hgvs.intervalmapper.Interval(s2, e2)
        intervalpair = hgvs.intervalmapper.IntervalPair(iv1, iv2)
        repr_str = 'IntervalPair(ref=Interval(start_i={},end_i={}),tgt=Interval(start_i={},end_i={}))'
        self.assertEqual(intervalpair.__repr__(), repr_str.format(s1, e1, s2, e2))

    @staticmethod
    def _build_mock_intervalmapper():
        iv1 = [hgvs.intervalmapper.Interval(1, 5), hgvs.intervalmapper.Interval(5, 6),
               hgvs.intervalmapper.Interval(6, 10)]

        iv2 = [hgvs.intervalmapper.Interval(1, 5), hgvs.intervalmapper.Interval(5, 5),
               hgvs.intervalmapper.Interval(5, 9)]

        ivp = [hgvs.intervalmapper.IntervalPair(iv1[i], iv2[i]) for i in xrange(len(iv1))]
        ivm = hgvs.intervalmapper.IntervalMapper(ivp)
        return ivm


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
