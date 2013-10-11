# -*- encoding: utf-8 -*-

import recordtype

class HGVSPosition( recordtype.recordtype('HGVSPosition', ['seqref','type','pos']) ):
    """
    HGVSPosition -- Represent partial HGVS tags that refer to a position without alleles
    """
    def __str__(self):
        return '{self.seqref}:{self.type}.{self.pos}'.format(self=self)
