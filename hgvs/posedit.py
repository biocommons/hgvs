import recordtype

from hgvs.exceptions import HGVSError

class PosEdit( recordtype.recordtype( 'PosEdit', [('pos',None),('edit',None),('uncertain',False)] )):
    """
    represents a **simple** variant
    """
    
    def __str__(self):
        rv = str(self.edit) if self.pos is None else '{self.pos}{self.edit}'.format(self=self)
        if self.uncertain:
            if self.edit == '?':
                pass
            elif self.edit == '0':
                rv = rv + '?'
            else:
                rv = '(' + rv + ')'
        return rv

    def set_uncertain(self):
        self.uncertain = True
        return self
    
