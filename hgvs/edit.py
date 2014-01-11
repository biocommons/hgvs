__doc__ = """
hgvs.edit -- representation of edit operations in HGVS variants

NARefAlt and AARefAlt are abstractions of several major variant
types.  They are distinguished by whether the ref and alt elements
of the structure.  The HGVS grammar for NA and AA are subtly
different (e.g., the ref AA in a protein substitution is part of the
location).


  TYPE      NA              AA          REF/ALT representation (NA; AA)
  subst     A>T             T           A/T           ''/T
  delins    del(AA)insTT    delinsTT    (AA|'')/TT    same
  del       del(AA)         del(AA)     (AA|'')/None  same
  ins       insTT           insTT       None/TT       same

patterns of seq/seq, seq/None, None/seq determine basic variant type
'' means missing optional sequence

"""

import recordtype

from hgvs.exceptions import HGVSError
from hgvs.utils import aa_to_aa1,aa1_to_aa3

class Edit(object):
    pass

class NARefAlt( Edit, recordtype.recordtype('NARefAlt', [('ref',None),('alt',None),('uncertain',False)]) ):

    def __str__(self):
        if self.ref is None and self.alt is None:
            raise HGVSError('RefAlt: ref and alt sequences are both undefined')

        # subst and delins
        if self.ref is not None and self.alt is not None:
            if self.ref == self.alt:
                s =  '='
            elif len(self.alt) == 1 and len(self.ref) == 1:
                s = '{self.ref}>{self.alt}'.format(self=self)
            else:
                s = 'del{self.ref}ins{self.alt}'.format(self=self)

        # del case
        elif self.ref is not None:
            s = 'del{self.ref}'.format(self=self)

        # ins case
        else:   # self.alt is not None
            s = 'ins{self.alt}'.format(self=self)

        return '('+s+')' if self.uncertain else s

    def set_uncertain(self):
        self.uncertain = True
        return self

        
class AARefAlt( Edit, recordtype.recordtype('AARefAlt', [('ref',None),('alt',None),('fs',None),('uncertain',False)]) ):
    def __init__(self,ref,alt,fs=None,uncertain=False):
        super(AARefAlt,self).__init__(ref=aa_to_aa1(ref),alt=aa_to_aa1(alt),fs=fs)

    def __str__(self):
        if self.ref is None and self.alt is None:
            #raise HGVSError('RefAlt: ref and alt sequences are both undefined')
            return '='

        # subst and delins
        if self.ref is not None and self.alt is not None:
            if self.ref == self.alt:
                s = '='
            elif len(self.ref) == 1 and len(self.alt) == 1:
                s = aa1_to_aa3(self.alt) + (self.fs or '')
            else:
                s = 'delins{alt}{fs}'.format(alt = aa1_to_aa3(self.alt),fs = self.fs or '')

        # del case
        elif self.ref is not None and self.alt is None:
            s = 'del'

        # ins case
        elif self.ref is None and self.alt is not None:
            s = 'ins{alt}{fs}'.format(alt=aa1_to_aa3(self.alt),fs=self.fs or '')

        else:
            raise RuntimeError("Should not be here")
            
        return '('+s+')' if self.uncertain else s
    
    def set_uncertain(self):
        self.uncertain = True
        return self


class AASub( AARefAlt ):
    def __str__(self):
        s = aa1_to_aa3(self.alt) + (self.fs or '')
        return '('+s+')' if self.uncertain else s


class AASpecial( Edit, recordtype.recordtype('AASpecial', [('status',None),('uncertain',False)]) ):
    def __str__(self):
        return '('+self.status+')' if self.uncertain else self.status
    
    def set_uncertain(self):
        self.uncertain = True
        return self


class Dup( Edit, recordtype.recordtype('Dup', [('seq',None),('uncertain',False)]) ):

    def __str__(self):
        return 'dup' + (self.seq or '')

    def set_uncertain(self):
        self.uncertain = True
        return self


class Repeat( Edit, recordtype.recordtype('Repeat', [('seq',None),('min',None),('max',None),('uncertain',False)]) ):

    def __str__(self):
        if self.min > self.max:
            raise HGVSError('Repeat min count must be less than or equal to max count')
        if self.min == self.max:
            return '{self.seq}[{self.min}]'.format(self=self)
        return '{self.seq}({self.min}_{self.max})'.format(self=self)

    def set_uncertain(self):
        self.uncertain = True
        return self





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


## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/invitae/hgvs)
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
