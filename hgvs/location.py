import recordtype


class Position(int):
    def __init__(self,position):
        super(Position,self).__init__()
        assert position>=1, self.__class__.__name__ + ': position must be >= 1'

class CDSPosition( recordtype.recordtype(
        'CDSPosition', field_names = [ 'base', 'offset' ] ) ):
    """
    Class for dealing with CDS coordinates in transcript variants.  Positions and
    offsets are 1-based, with no 0, per the HGVS recommendations.  (If you're using this with UTA, 
    be aware that UTA uses interbase coordinates.)

    This class models CDS positions using a `base' coordinate and an `offset' with
    the following interpretations:

    position   base     offset     meaning
    c.55        55         0       cds position 55
    c.55+1      55         1       intronic variant +1 from boundary
    c.-55        0       -55       5' UTR variant, 55 nt upstream of ATG
    c.*55        0        55       3' UTR variant, 55 nt after STOP
    c.1          1         0	   start codon
    c.1234    1234         0       stop codon (assuming CDS length is 1234)
    c.*1         0         1       STOP + 1

    In other words:
    - offset == 0 for coding positions, and offset != 0 for non-coding positions.
    - base == 0 for UTR, with offset<0 for 5' and offset>0 for 3'
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


class Interval( recordtype.recordtype(
        'Interval', field_names = [ 'start', 'end' ] ) ):
    def __init__(self,start,end):
        super(Interval,self).__init__(start,end)

    def __str__(self):
        if self.end is None or self.start == self.end:
            return str(self.start)
        return str(self.start) + '_' + str(self.end)

class CDSInterval( Interval ):
    def __init__(self,start,end):
        super(CDSInterval,self).__init__(start,end)
        if end is not None:
            assert self.end.base >= self.start.base, 'end must be >= than start'
            assert not (self.start.base == self.end.base and self.start.offset >= self.end.offset), 'end offset must be >= than start offset'
    
