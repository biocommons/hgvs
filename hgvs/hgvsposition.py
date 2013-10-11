# -*- encoding: utf-8 -*-
__doc__ = """
"""

import recordtype

class HGVSPosition( recordtype.recordtype('HGVSPosition', ['seqref','type','pos']) ):
    """
    """
    
    def __str__(self):
        return '{self.seqref}:{self.type}.{self.posedits}'.format(self=self)
