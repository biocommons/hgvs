# -*- encoding: utf-8 -*-
__doc__ = """
hgvs.location -- classes for dealing with the locations of HGVS variants

This module provides for Representing the location of variants in HGVS nomenclature, including:

  * integers and integer intervals (e.g., NC_012345.6:g.3403243_3403248A>C)
  * CDS positions and intervals (e.g., NM_01234.5:c.56+12_56+14delAC)
  * CDS stop coordinates (e.g., NM_01234.5:c.*13A>C)  

Classes:
  * SimplePosition -- a simple integer
  * BaseOffsetPosition -- a position with datum, base, and offset for c. and r. coordinates
  * AAPosition -- an amino acid position (with AA)
  * Interval -- an interval of Positions
"""

import recordtype

from hgvs.utils import aa1_to_aa3,aa_to_aa1

SEQ_START = 0
CDS_START = 1
CDS_END = 2

class SimplePosition(int):
    def __init__(self,position):
        super(SimplePosition,self).__init__()
        assert position>=1, self.__class__.__name__ + ': position must be >= 1'

    @property
    def base(self):
        """convenience property for similarity with BaseOffsetPosition"""
        return self


class BaseOffsetPosition( recordtype.recordtype(
        'BaseOffsetPosition', field_names = [ 'base', 'offset', 'datum' ] ) ):
    """
    Class for dealing with CDS coordinates in transcript variants.

    This class models CDS positions using a `base' coordinate, which is
    measured relative to a specified `datum' (CDS_START or CDS_END), and
    an `offset', which is 0 for exonic positions and non-zero for intronic
    positions.  Positions and offsets are 1-based, with no 0, per the HGVS
    recommendations.  (If you're using this with UTA, be aware that UTA
    uses interbase coordinates.)

    hgvs     datum      base  offset  meaning
    r.55     SEQ_START    55       0  RNA position 55
    c.55     CDS_START    55       0  CDS position 55
    c.55     CDS_START    55       0  CDS position 55
    c.55+1   CDS_START    55       1  intronic variant +1 from boundary
    c.-55    CDS_START   -55       0  5' UTR variant, 55 nt upstream of ATG
    c.1      CDS_START     1       0    start codon
    c.1234   CDS_START  1234       0  stop codon (assuming CDS length is 1233)
    c.*1     CDS_END       0       1  STOP + 1
    c.*55    CDS_END       0      55  3' UTR variant, 55 nt after STOP
    """
    
    def __init__(self,base,offset=0,datum=SEQ_START):
        assert base != 0, 'BaseOffsetPosition base may not be 0'
        assert datum == CDS_START or base >= 1, 'BaseOffsetPosition base must be >=1 for datum = SEQ_START or CDS_END'
        super(BaseOffsetPosition,self).__init__(base=base,offset=offset,datum=datum)

    def __str__(self):
        base_str = '*' + str(self.base) if self.datum == CDS_END else str(self.base)
        offset_str = '' if self.offset == 0 else '%+d' % self.offset
        return base_str + offset_str


class AAPosition( recordtype.recordtype(
        'AAPosition', field_names = [ 'pos', 'aa' ] ) ):
    def __init__(self,pos,aa):
        super(AAPosition,self).__init__(pos,aa_to_aa1(aa))

    def __str__(self):
        return aa1_to_aa3(self.aa) + str(self.pos)


class Interval( recordtype.recordtype(
        'Interval', field_names = [ 'start', 'end' ] ) ):
    def __init__(self,start,end):
        super(Interval,self).__init__(start,end)

    def __str__(self):
        if self.end is None or self.start == self.end:
            return str(self.start)
        return str(self.start) + '_' + str(self.end)
