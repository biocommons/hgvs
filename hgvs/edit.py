__doc__ = """
hgvs.edit -- representation of edit operations in HGVS variants

"""

import recordtype

from hgvs.exceptions import HGVSError

class Edit(object):
    pass

class DelIns( Edit, recordtype.recordtype('DelIns', ['ref','alt'], default=None) ):
    """
    DelIns is an abstraction of several major variant types.  They are distinguished by 
    whether the ref and alt elements of the structure.

    TYPE                                        REF             ALT
    SNV (e.g., c.123A>T)                        single nt       single nt
    Deletion (c.123del)                         sequence        None
    Insertion (c.123_124insT)                   None            sequence
    DeletionInsertion (c.123_130del8insAT)      sequence        sequence
    """
    def __str__(self):
        if self.ref is None and self.alt is None:
            raise HGVSError('DelIns: ref and alt sequences are both empty')
        if self.ref is not None and self.alt is not None:
            if self.ref == self.alt:
                return '='
            if len(self.ref) == 1 and len(self.alt) == 1:
                return '{self.ref}>{self.alt}'.format(self=self)
            return 'del{self.ref}ins{self.alt}'.format(self=self)
        if self.ref is not None and self.alt is None:
            return 'del{self.ref}'.format(self=self)
        if self.ref is None and self.alt is not None:
            return 'ins{self.alt}'.format(self=self)

class Dup( Edit, recordtype.recordtype('Dup', ['seq'], default=None) ):
    def __str__(self):
        return 'dup' + (self.seq or '')

class Repeat( Edit, recordtype.recordtype('Repeat', ['seq','min','max'], default=None) ):
    def __str__(self):
        if self.min > self.max:
            raise HGVSError('Repeat min count must be less than or equal to max count')
        if self.min == self.max:
            return '{self.seq}[{self.min}]'.format(self=self)
        return '{self.seq}({self.min}_{self.max})'.format(self=self)

# class Inv( Edit, recordtype.recordtype('Inv', [], default=None) ):
#     def __str__(self):
#         return ''
# 
# class Con( Edit, recordtype.recordtype('Con', ['con'], default=None) ):
#     def __str__(self):
#         return self.con
# 
# class ComplexVariant( Edit, recordtype.recordtype('ComplexVariant', ['edits','rel'], default=None) ):
#     def __str__(self):
#         return '[' + self.rel.join( self.edits ) + ']'
# 
# class CompoundVariant( Edit, recordtype.recordtype('CompoundVariant', ['edits'], default=None) ):
#     def __str__(self):
#         return ';'.join( [ '['+e+']' for e in self.edits ] )
# 
# class MosaicVariant( Edit, recordtype.recordtype('Edit', ['edit'], default=None) ):
#     def __str__(self):
#         return '[=/{self.edit}]'.format(self=self)
# 
# class ChimericVariant( Edit, recordtype.recordtype('Edit', ['edit'], default=None) ):
#     def __str__(self):
#         return '[=//{self.edit}]'.format(self=self)

