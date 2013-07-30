import recordtype


class Position(int):
    def __init__(self,position):
        super(Position,self).__init__()
        assert position>=1, self.__class__.__name__ + ': position must be >= 1'

class CDSPosition( recordtype.recordtype(
        'CDSPosition', field_names = [ 'base', 'offset' ] ) ):
    """
    Class for dealing with CDS coordinates in transcript variants.  Positions and
    offsets are 1-based, per the HGVS recommendations.

    This class models CDS positions using a `base' coordinate and an `offset' with
    the following interpretations:

    position   base     offset     meaning
    c.55        55         0       cds position 55
    c.55+1      55         1       intronic variant +1 from boundary
    c.-55        0       -55       5' UTR variant
    c.*55        0        55       3' UTR variant (* is STOP)

    In other words:
    - offset == 0 for coding positions, and offset != 0 for non-coding positions.
    - base == 0 for UTR
    """
    
    def __init__(self,base,offset=0):
        super(CDSPosition,self).__init__(base=base,offset=offset)

    @property
    def is_coding(self):
        return self.offset == 0

    @property
    def is_utr(self):
        return self.base == 0

    @property
    def is_intronic(self):
        return self.base != 0 and self.offset != 0

    def __str__(self):
        if self.offset == 0:              # coding
            return str(self.base)
        elif self.base != 0:              # intronic
            return '%d%+d' % (self.base, self.offset)
        elif self.offset < 0:             # 5' UTR
            return str(self.offset)
        elif self.offset > 0:             # 3' UTR
            return '*' + str(self.offset)
        raise RuntimeError('Fell through to unexpected <base,offset> case')


class Range( recordtype.recordtype(
        'Range', field_names = [ 'start', 'end' ] ) ):
    def __init__(self,start,end):
        super(Range,self).__init__(start,end)
        if end is not None:
            assert self.end.base >= self.start.base, 'end must be >= than start'
            assert not (self.start.base == self.end.base and self.start.offset >= self.end.offset), 'end offset must be >= than start offset'

    def __str__(self):
        if self.end is None or self.start == self.end:
            return str(self.start)
        return str(self.start) + '_' + str(self.end)
