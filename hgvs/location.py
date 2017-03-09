# -*- coding: utf-8 -*-
"""Provides classes for dealing with the locations of HGVS variants

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

from __future__ import absolute_import, division, print_function, unicode_literals

import recordtype

from functools import total_ordering
from bioutils.sequences import aa1_to_aa3

import hgvs
from hgvs.exceptions import HGVSUnsupportedOperationError, HGVSInvalidIntervalError
from hgvs.enums import ValidationLevel


SEQ_START = 0
CDS_START = 1
CDS_END = 2


@total_ordering
class SimplePosition(recordtype.recordtype("SimplePosition", field_names=[("base", None), ("uncertain", False)])):
    def __str__(self):
        self.validate()
        s = "?" if self.base is None else str(self.base)
        return "(" + s + ")" if self.uncertain else s

    def format(self, conf):
        return str(self)

    @property
    def is_uncertain(self):
        """return True if the position is marked uncertain or undefined"""
        return self.uncertain or self.base in None

    def _set_uncertain(self):
        "mark this location as uncertain and return reference to self; this is called during parsing (see hgvs.ometa)"
        self.uncertain = True
        return self

    def validate(self):
        if self.base is not None and self.base < 1:
            return (ValidationLevel.ERROR, "Position base must be >= 1")
        return (ValidationLevel.VALID, None)

    def __sub__(lhs, rhs):
        assert type(lhs) == type(rhs), "Cannot substract coordinates of different representations"
        return lhs.base - rhs.base

    def __eq__(lhs, rhs):
        assert type(lhs) == type(rhs), "Cannot compare coordinates of different representations"
        if lhs.uncertain or rhs.uncertain:
            raise HGVSUnsupportedOperationError("Cannot compare coordinates of uncertain positions")
        return lhs.base == rhs.base

    def __lt__(lhs, rhs):
        assert type(lhs) == type(rhs), "Cannot compare coordinates of different representations"
        if lhs.uncertain or rhs.uncertain:
            raise HGVSUnsupportedOperationError("Cannot compare coordinates of uncertain positions")
        return lhs.base < rhs.base


@total_ordering
class BaseOffsetPosition(
        recordtype.recordtype(
            'BaseOffsetPosition',
            field_names=[('base', None), ('offset', 0), ('datum', SEQ_START), ('uncertain', False)])):
    """
    Class for dealing with CDS coordinates in transcript variants.

    This class models CDS positions using a `base` coordinate, which is
    measured relative to a specified `datum` (CDS_START or CDS_END), and
    an `offset`, which is 0 for exonic positions and non-zero for intronic
    positions.  **Positions and offsets are 1-based**, with no 0, per the HGVS
    recommendations.  (If you"re using this with UTA, be aware that UTA
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
    | c.1      | CDS_START  |    1  |      0  | start codon                              |
    +----------+------------+-------+---------+------------------------------------------+
    | c.1234   | CDS_START  | 1234  |      0  | stop codon (assuming CDS length is 1233) |
    +----------+------------+-------+---------+------------------------------------------+
    | c.*1     | CDS_END    |    0  |      1  | STOP + 1                                 |
    +----------+------------+-------+---------+------------------------------------------+
    | c.*55    | CDS_END    |    0  |     55  | 3' UTR variant, 55 nt after STOP         |
    +----------+------------+-------+---------+------------------------------------------+
    """

    def validate(self):
        if self.base is not None and self.base == 0:
            return (ValidationLevel.ERROR, "BaseOffsetPosition base may not be 0")
        if self.base is not None and self.datum != CDS_START and self.base < 1:
            return (ValidationLevel.ERROR, "BaseOffsetPosition base must be >=1 for datum = SEQ_START or CDS_END")
        return (ValidationLevel.VALID, None)

    def __str__(self):
        self.validate()
        base_str = ("?" if self.base is None else "*" + str(self.base) if self.datum == CDS_END else str(self.base))
        offset_str = ("+?" if self.offset is None else "" if self.offset == 0 else "%+d" % self.offset)
        pos = base_str + offset_str
        return "(" + pos + ")" if self.uncertain else pos

    def format(self, conf):
        return str(self)

    def _set_uncertain(self):
        "mark this location as uncertain and return reference to self; this is called during parsing (see hgvs.ometa)"
        self.uncertain = True
        return self

    @property
    def is_uncertain(self):
        """return True if the position is marked uncertain or undefined"""
        return self.uncertain or self.base is None or self.offset is None

    @property
    def is_intronic(self):
        """returns True if the variant is intronic (if the offset is None or non-zero)"""
        return (self.offset is None or self.offset != 0)

    def __sub__(lhs, rhs):
        assert type(lhs) == type(rhs), "Cannot substract coordinates of different representations"
        if lhs.datum != rhs.datum:
            raise HGVSUnsupportedOperationError("Interval length measured from different datums is ill-defined")
        if lhs.base == rhs.base:
            return lhs.offset - rhs.offset
        if lhs.offset != 0 or rhs.offset != 0:
            raise HGVSUnsupportedOperationError("Interval length with intronic offsets is ill-defined")
        straddles_zero = 1 if (lhs.base > 0 and rhs.base < 0) else 0
        return lhs.base - rhs.base - straddles_zero

    def __eq__(lhs, rhs):
        assert type(lhs) == type(rhs), "Cannot compare coordinates of different representations"
        if lhs.uncertain or rhs.uncertain:
            raise HGVSUnsupportedOperationError("Cannot compare coordinates of uncertain positions")
        return lhs.datum == rhs.datum and lhs.base == rhs.base and lhs.offset == rhs.offset

    def __lt__(lhs, rhs):
        assert type(lhs) == type(rhs), "Cannot compare coordinates of different representations"
        if lhs.uncertain or rhs.uncertain:
            raise HGVSUnsupportedOperationError("Cannot compare coordinates of uncertain positions")
        if lhs.datum == rhs.datum:
            if lhs.base == rhs.base:
                return lhs.offset < rhs.offset
            else:
                if ((rhs.base - lhs.base == 1 and lhs.offset > 0 and rhs.offset < 0) or
                    (lhs.base - rhs.base == 1 and rhs.offset > 0 and lhs.offset < 0)):
                    raise HGVSUnsupportedOperationError(
                        "Cannot compare coordinates in the same intron with one based on end of exon and the other based on start of next exon"
                    )
                else:
                    return lhs.base < rhs.base
        else:
            if lhs.datum == SEQ_START or rhs.datum == SEQ_START:
                raise HGVSUnsupportedOperationError(
                    "Cannot compare coordinates of datum SEQ_START with CDS_START or CDS_END")
            else:
                return lhs.datum < rhs.datum


class AAPosition(recordtype.recordtype("AAPosition", field_names=[("base", None), ("aa", None), ("uncertain", False)])):
    def validate(self):
        if self.base is not None and self.base < 1:
            return (ValidationLevel.ERROR, "AAPosition location must be >=1")
        if len(self.aa) != 1:
            return (ValidationLevel.ERROR, "More than 1 AA associated with position")
        return (ValidationLevel.VALID, None)

    def format(self, conf=None):
        self.validate()

        p_3_letter = hgvs.global_config.formatting.p_3_letter
        p_term_asterisk = hgvs.global_config.formatting.p_term_asterisk
        if conf and "p_3_letter" in conf and conf["p_3_letter"] is not None:
            p_3_letter = conf["p_3_letter"]
        if conf and "p_term_asterisk" in conf and conf["p_term_asterisk"] is not None:
            p_term_asterisk = conf["p_term_asterisk"]

        pos = "?" if self.base is None else str(self.base)
        if p_3_letter:
            aa = "?" if self.aa is None else aa1_to_aa3(self.aa)
            if p_term_asterisk and aa == "Ter":
                aa = "*"
        else:
            aa = "?" if self.aa is None else self.aa
        s = aa + pos
        return "(" + s + ")" if self.uncertain else s

    __str__ = format

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

    def __sub__(lhs, rhs):
        assert type(lhs) == type(rhs), "Cannot substract coordinates of different representations"
        return lhs.base - rhs.base

    def __lt__(lhs, rhs):
        assert type(lhs) == type(rhs), "Cannot compare coordinates of different representations"
        if lhs.uncertain or rhs.uncertain:
            raise HGVSUnsupportedOperationError("Cannot compare coordinates of uncertain positions")
        return lhs.base < rhs.base

    def __gt__(lhs, rhs):
        assert type(lhs) == type(rhs), "Cannot compare coordinates of different representations"
        if lhs.uncertain or rhs.uncertain:
            raise HGVSUnsupportedOperationError("Cannot compare coordinates of uncertain positions")
        return lhs.base > rhs.base


class Interval(recordtype.recordtype("Interval", field_names=["start", ("end", None), ("uncertain", False)])):
    def validate(self):
        if self.start:
            (res, msg) = self.start.validate()
            if res != ValidationLevel.VALID:
                return (res, msg)
        if self.end:
            (res, msg) = self.end.validate()
            if res != ValidationLevel.VALID:
                return (res, msg)
        # Check start less than or equal to end
        if not self.start or not self.end:
            return (ValidationLevel.VALID, None)
        try:
            if self.start <= self.end:
                return (ValidationLevel.VALID, None)
            else:
                return (ValidationLevel.ERROR, "base start position must be <= end position")
        except HGVSUnsupportedOperationError as err:
            return (ValidationLevel.WARNING, str(err))

    def format(self, conf=None):
        if self.end is None or self.start == self.end:
            return self.start.format(conf)
        iv = self.start.format(conf) + "_" + self.end.format(conf)
        return "(" + iv + ")" if self.uncertain else iv

    __str__ = format

    def _set_uncertain(self):
        "mark this interval as uncertain and return reference to self; this is called during parsing (see hgvs.ometa)"
        self.uncertain = True
        return self

    def _length(self):
        return 1 if self.end is None else self.end - self.start + 1

    @property
    def is_uncertain(self):
        """return True if the position is marked uncertain or undefined"""
        return self.uncertain or self.start.is_uncertain or self.end.is_uncertain


class BaseOffsetInterval(Interval):
    """BaseOffsetInterval isa Interval of BaseOffsetPositions.  The only
    additional functionality over Interval is to ensure that the dutum
    of end and start are compatible.

    """

    def __init__(self, *args, **kwargs):
        super(BaseOffsetInterval, self).__init__(*args, **kwargs)

        # #330: In a post-ter interval like *87_91, the * binds only
        # to the start. This means that the start.datum is CDS_END,
        # but the end.datum is CDS_START (the default). 
        if self.start.datum == CDS_END:
            self.end.datum = CDS_END

        self.check_datum()

    def check_datum(self):
        # check for valid combinations of start and end datums
        if (self.start.datum, self.end.datum) not in [
            (SEQ_START, SEQ_START),
            (CDS_START, CDS_START),
            (CDS_START, CDS_END),
            (CDS_END, CDS_END),
        ]:
            raise HGVSInvalidIntervalError("BaseOffsetInterval start datum and end datum are incompatible")


# <LICENSE>
# Copyright 2013-2015 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
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
