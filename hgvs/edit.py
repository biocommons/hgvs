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
                s = '='
            elif len(self.alt) == 1 and len(self.ref) == 1 and not self.ref.isdigit(): # don't turn del5insT into 5>T
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

    @property
    def type(self):
        if self.ref is not None and self.alt is not None:
            if self.ref == self.alt:
                edit_type = None
            elif len(self.alt) == 1 and len(self.ref) == 1 and not self.ref.isdigit():
                edit_type = 'sub'
            else:
                edit_type = 'delins'
        elif self.ref is not None:
            edit_type = 'del'
        else:
            edit_type = 'ins'
        return edit_type

class AARefAlt( Edit, recordtype.recordtype('AARefAlt', [('ref',None),('alt',None), ('uncertain',False)]) ):
    def __init__(self,ref, alt, uncertain=False):
        super(AARefAlt, self).__init__(ref=aa_to_aa1(ref), alt=aa_to_aa1(alt), uncertain=uncertain)

    def __str__(self):
        if self.ref is None and self.alt is None:
            #raise HGVSError('RefAlt: ref and alt sequences are both undefined')
            return '='

        # subst and delins
        if self.ref is not None and self.alt is not None:
            if self.ref == self.alt:
                s = '='
            elif len(self.ref) == 1 and len(self.alt) == 1:
                s = aa1_to_aa3(self.alt)
            else:
                s = 'delins{alt}'.format(alt = aa1_to_aa3(self.alt))

        # del case
        elif self.ref is not None and self.alt is None:
            s = 'del'

        # ins case
        elif self.ref is None and self.alt is not None:
            s = 'ins{alt}'.format(alt=aa1_to_aa3(self.alt))

        else:
            raise RuntimeError("Should not be here")
            
        return '('+s+')' if self.uncertain else s
    
    def set_uncertain(self):
        self.uncertain = True
        return self

    @property
    def type(self):
        if self.ref is not None and self.alt is not None:
            if self.ref == self.alt:
                edit_type = None
            elif len(self.ref) == 1 and len(self.alt) == 1:
                edit_type = 'sub'
            else:
                edit_type = 'delins'
        elif self.ref is not None and self.alt is None:
            edit_type = 'del'
        elif self.ref is None and self.alt is not None:
            edit_type = 'ins'
        return edit_type

class AASub( AARefAlt ):
    def __str__(self):
        s = aa1_to_aa3(self.alt) if self.alt != '?' else self.alt
        return '('+s+')' if self.uncertain else s

    @property
    def type(self):
        return 'sub'

class AAFs(Edit, recordtype.recordtype('AAFs', [('ref',None),('alt',None),('length',None),('uncertain',False)])):
    def __init__(self,ref,alt,length=None,uncertain=False):
        super(AAFs, self).__init__(ref=aa_to_aa1(ref), alt=aa_to_aa1(alt), length=length, uncertain=uncertain)

    def __str__(self):
        st_length = self.length or ''
        s = "{alt}fsTer{length}".format(alt=aa1_to_aa3(self.alt), length=st_length)
        return '('+s+')' if self.uncertain else s

    def set_uncertain(self):
        self.uncertain = True
        return self

    @property
    def type(self):
        return 'fs'

class AAExt(Edit, recordtype.recordtype('AAExt', [('ref',None),('alt',None), ('aaterm', None), ('length',None),
                                                  ('uncertain',False)])):
    def __init__(self,ref,alt,aaterm=None, length=None,uncertain=False):
        super(AAExt, self).__init__(ref=aa_to_aa1(ref), alt=aa_to_aa1(alt), aaterm=aa_to_aa1(aaterm), length=length,
                                    uncertain=uncertain)

    def __str__(self):
        st_alt = self.alt or ''
        st_aaterm = self.aaterm or ''
        st_length = self.length or ''
        s = "{alt}ext{term}{length}".format(alt=aa1_to_aa3(st_alt), term=aa1_to_aa3(st_aaterm), length=st_length)
        return '('+s+')' if self.uncertain else s

    def set_uncertain(self):
        self.uncertain = True
        return self

    @property
    def type(self):
        return 'ext'

class Dup( Edit, recordtype.recordtype('Dup', [('seq',None),('uncertain',False)]) ):

    def __str__(self):
        return 'dup' + (self.seq or '')

    def set_uncertain(self):
        self.uncertain = True
        return self

    @property
    def type(self):
        return 'dup'

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

    @property
    def type(self):
        return 'repeat'




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
