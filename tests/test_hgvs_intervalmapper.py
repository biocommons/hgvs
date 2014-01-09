from __future__ import with_statement

import unittest

import hgvs.exceptions
import hgvs.intervalmapper

class Test_IntervalMapper(unittest.TestCase):
    longMessage = True

    def setUp(self):
        pass

    def test_intervalmapper_valid(self):
        interval = hgvs.intervalmapper.Interval(3, 4)
        st = str(interval)
        self.assertEqual(1, interval.len)
        self.assertEqual(interval.__repr__(), 'Interval(start_i=3,end_i=4)')

    def test_intervalmapper_invalid_start_gt_end(self):
        with self.assertRaises(hgvs.exceptions.InvalidIntervalError):
            interval = hgvs.intervalmapper.Interval(4,3)

    def test_intervalpair_valid_match(self):
        iv1 = hgvs.intervalmapper.Interval(3, 4)
        iv2 = hgvs.intervalmapper.Interval(6, 7)
        intervalpair = hgvs.intervalmapper.IntervalPair(iv1, iv2)
        self.assertEqual(intervalpair.__repr__(),
                         'IntervalPair(ref=Interval(start_i=3,end_i=4),tgt=Interval(start_i=6,end_i=7))')

    def test_intervalpair_valid_insertion(self):
        iv1 = hgvs.intervalmapper.Interval(3, 3)
        iv2 = hgvs.intervalmapper.Interval(6, 7)
        intervalpair = hgvs.intervalmapper.IntervalPair(iv1, iv2)
        self.assertEqual(intervalpair.__repr__(),
                         'IntervalPair(ref=Interval(start_i=3,end_i=3),tgt=Interval(start_i=6,end_i=7))')

    def test_intervalpair_valid_deletion(self):
        iv1 = hgvs.intervalmapper.Interval(3, 4)
        iv2 = hgvs.intervalmapper.Interval(6, 6)
        intervalpair = hgvs.intervalmapper.IntervalPair(iv1, iv2)
        self.assertEqual(intervalpair.__repr__(),
                         'IntervalPair(ref=Interval(start_i=3,end_i=4),tgt=Interval(start_i=6,end_i=6))')

    def test_intervalpair_invalid(self):
        iv1 = hgvs.intervalmapper.Interval(3, 4)
        iv2 = hgvs.intervalmapper.Interval(6, 9)
        with self.assertRaises(hgvs.exceptions.InvalidIntervalError):
            intervalpair = hgvs.intervalmapper.IntervalPair(iv1, iv2)


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
