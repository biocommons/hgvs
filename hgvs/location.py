# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
"""
hgvs.location -- classes for dealing with the locations of HGVS variants

This module provides for Representing the location of variants in HGVS nomenclature, including:

  * integers and integer intervals (e.g., NC_012345.6:g.3403243_3403248A>C)
  * CDS positions and intervals (e.g., NM_01234.5:c.56+12_56+14delAC)
  * CDS stop coordinates (e.g., NM_01234.5:c.*13A>C)  

Classes:

  * :class:`SimplePosition` -- a simple integer
  * :class:`BaseOffsetPosition` -- a position with datum, base, and offset for c. and r. coordinates
  * :class:`AAPosition` -- an amino acid position (with AA)
  * :class:`Interval` -- an interval of Positions
"""

import recordtype

from bioutils.sequences import aa1_to_aa3

SEQ_START = 0
CDS_START = 1
CDS_END = 2


class SimplePosition(recordtype.recordtype('SimplePosition', field_names=[('base', None), ('uncertain', False)])):
    def __str__(self):
        self.validate()
        s = '?' if self.base is None else str(self.base)
        return '(' + s + ')' if self.uncertain else s

    @property
    def is_uncertain(self):
        """return True if the position is marked uncertain or undefined"""
        return self.uncertain or self.base in None

    def _set_uncertain(self):
        "mark this location as uncertain and return reference to self; this is called during parsing (see hgvs.ometa)"
        self.uncertain = True
        return self

    def validate(self):
        "raise AssertionError if instance variables are invalid; otherwise return True"
        assert self.base is None or self.base >= 1, self.__class__.__name__ + ': base must be >= 1'
        return True


class BaseOffsetPosition(recordtype.recordtype(
    'BaseOffsetPosition',
    field_names=[('base', None), ('offset', 0), ('datum', SEQ_START), ('uncertain', False)])):
    """
    Class for dealing with CDS coordinates in transcript variants.

    This class models CDS positions using a `base` coordinate, which is
    measured relative to a specified `datum` (CDS_START or CDS_END), and
    an `offset`, which is 0 for exonic positions and non-zero for intronic
    positions.  **Positions and offsets are 1-based**, with no 0, per the HGVS
    recommendations.  (If you're using this with UTA, be aware that UTA
    uses interbase coordinates.)

    +----------+------------+-------+---------+------------------------------------------+
    | hgvs     | datum      | base  | offset  | meaning                                  |
    +==========+============+=======+=========+==========================================+
    | r.55     | SEQ_START  |   55  |      0  | RNA position 55                          |
    +----------+------------+-------+---------+------------------------------------------+
    | c.55     | CDS_START  |   55  |      0  | CDS position 55                          |
    +----------+------------+-------+---------+------------------------------------------+
    | c.55     | CDS_START  |   55  |      0  | CDS position 55                          |
    +----------+------------+-------+---------+------------------------------------------+
    | c.55+1   | CDS_START  |   55  |      1  | intronic variant +1 from boundary        |
    +----------+------------+-------+---------+------------------------------------------+
    | c.-55    | CDS_START  |  -55  |      0  | 5' UTR variant, 55 nt upstream of ATG    |
    +----------+------------+-------+---------+------------------------------------------+
    | c.1      | CDS_START  |    1  |      0  |   start codon                            |
    +----------+------------+-------+---------+------------------------------------------+
    | c.1234   | CDS_START  | 1234  |      0  | stop codon (assuming CDS length is 1233) |
    +----------+------------+-------+---------+------------------------------------------+
    | c.*1     | CDS_END    |    0  |      1  | STOP + 1                                 |
    +----------+------------+-------+---------+------------------------------------------+
    | c.*55    | CDS_END    |    0  |     55  | 3' UTR variant, 55 nt after STOP         |
    +----------+------------+-------+---------+------------------------------------------+
    """

    def validate(self):
        "raise AssertionError if instance variables are invalid; otherwise return True"
        assert self.base is None or self.base != 0, 'BaseOffsetPosition base may not be 0'
        assert self.base is None or self.datum == CDS_START or self.base >= 1, 'BaseOffsetPosition base must be >=1 for datum = SEQ_START or CDS_END'
        return True

    def __str__(self):
        self.validate()
        base_str = ('?' if self.base is None else '*' + str(self.base) if self.datum == CDS_END else str(self.base))
        offset_str = ('+?' if self.offset is None else '' if self.offset == 0 else '%+d' % self.offset)
        pos = base_str + offset_str
        return '(' + pos + ')' if self.uncertain else pos

    def _set_uncertain(self):
        "mark this location as uncertain and return reference to self; this is called during parsing (see hgvs.ometa)"
        self.uncertain = True
        return self

    @property
    def is_uncertain(self):
        """return True if the position is marked uncertain or undefined"""
        return self.uncertain or self.base is None or self.offset is None


class AAPosition(recordtype.recordtype('AAPosition', field_names=[('base', None), ('aa', None), ('uncertain', False)])):
    def validate(self):
        "raise AssertionError if instance variables are invalid; otherwise return True"
        assert self.base is None or self.base >= 1, 'AAPosition location must be >=1'
        assert len(self.aa) == 1, 'More than 1 AA associated with position'
        return True

    def __str__(self):
        self.validate()
        pos = '?' if self.base is None else str(self.base)
        aa = '?' if self.aa is None else aa1_to_aa3(self.aa)
        s = aa + pos
        return '(' + s + ')' if self.uncertain else s

    @property
    def pos(self):
        """return base, for backward compatibility"""
        return self.base

    def _set_uncertain(self):
        "mark this location as uncertain and return reference to self; this is called during parsing (see hgvs.ometa)"
        self.uncertain = True
        return self

    @property
    def is_uncertain(self):
        """return True if the position is marked uncertain or undefined"""
        return self.uncertain or self.base is None or self.aa is None


class Interval(recordtype.recordtype('Interval', field_names=['start', ('end', None), ('uncertain', False)])):
    def validate(self):
        "raise AssertionError if instance variables are invalid; otherwise return True"
        return True

    def __str__(self):
        self.validate()
        if self.end is None or self.start == self.end:
            return str(self.start)
        iv = str(self.start) + '_' + str(self.end)
        return '(' + iv + ')' if self.uncertain else iv

    def _set_uncertain(self):
        "mark this interval as uncertain and return reference to self; this is called during parsing (see hgvs.ometa)"
        self.uncertain = True
        return self

    @property
    def is_uncertain(self):
        """return True if the position is marked uncertain or undefined"""
        return self.uncertain or self.start.is_uncertain or self.end.is_uncertain

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
