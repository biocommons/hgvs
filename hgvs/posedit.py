import recordtype

from hgvs.exceptions import HGVSError

class PosEdit( recordtype.recordtype('PosEdit', ['pos','edit'], default=None) ):
    """
    represents a **simple** variant
    """
    
    def __str__(self):
        return '{self.pos}{self.edit}'.format(self=self)
