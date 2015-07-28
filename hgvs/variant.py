# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

from .decorators.deprecated import deprecated

from hgvs.exceptions import HGVSInvalidVariantError
import hgvs.posedit
import recordtype


class SequenceVariant(recordtype.recordtype('SequenceVariant', ['ac', 'type', 'posedit'])):
    """
    represents a basic HGVS variant.  The only requirement is that each
    component can be stringified; for example, passing pos as either a string
    or an hgvs.location.CDSInterval (for example) are both intended uses
    """

    def __str__(self):
        if self.ac is not None:
            return '{self.ac}:{self.type}.{self.posedit}'.format(self=self)
        else:
            return '{self.type}.{self.posedit}'.format(self=self)




PHASE1 = 0
PHASE2 = 1


class ComplexVariant(recordtype.recordtype('ComplexVariant', ['variants', 'phases', 'uncertain'])):
    """
    The abstract class for CompoundVariant, MosiacVariant and ChimericVariant,
    which are stored as a list of SequenceVariant and the phase info for each
    variant in the list
    """
    
    separator = ''
    
    def __init__(self, variants, phases, uncertain=False):
        super(ComplexVariant, self).__init__(variants, phases, uncertain)
        if not self.all_same_type:
            raise HGVSInvalidVariantError("ComplexVariant does not support variants "
                                          "with different types.")
        if self.phases and len(self.phases) != len(self.variants):
            raise HGVSInvalidVariantError("The number of edits and the number of phase "
                                          "info should be equal for ComplexVariant.")
    
    
    
    def __getitem__(self, i):
        return self.variants[i]

    def __setitem__(self, i, value):
        self.variants[i] = value
    
    def __iter__(self):
        for var in self.variants:
            yield var

    def __len__(self):
        return len(self.variants)
    
    def __str__(self):
        var_type = self.variants[0].type
        
        if self.uncertain:
            separator = '('+self.separator+')'
        else:
            separator = self.separator
        
        if self.phases:
            vars1 = []
            vars2 = []
            for var, phase in zip(self.variants, self.phases):
                if phase == PHASE1:
                    vars1.append(var)
                else:
                    vars2.append(var)
        
        if self.all_same_ac:
            var_ac = self.variants[0].ac
            
            if self.phases:
                vars1_posedits = separator.join([ str(var.posedit) for var in vars1 ])
                vars2_posedits = separator.join([ str(var.posedit) for var in vars2 ])
                var_posedits   = '[' + vars1_posedits + '];[' + vars2_posedits + ']'
            else:
                var_posedits = separator.join([ str(var.posedit) for var in self.variants ])
                if len(self.variants) > 1:
                    var_posedits = '[' + var_posedits + ']'
            
            if var_ac is not None:
                return '{ac}:{type}.{posedits}'.format(ac=var_ac, type=var_type, posedits=var_posedits)
            else:
                return '{type}.{posedits}'.format(type=var_type, posedits=var_posedits)
        else:
            if self.phases:
                vars1_var = separator.join(map(str, vars1))
                vars2_var = separator.join(map(str, vars2))
                vars      = '[' + vars1_var + '];[' + vars2_var + ']'
            else:
                vars = separator.join(map(str, self.variants))
                vars = '[' + vars + ']'
            
            return vars
    
    
    def _set_uncertain(self):
        self.uncertain = True
        return self
    
    @property
    def all_same_ac(self):
        """Return True if all variants have same accession or no accession, otherwise False."""
        num = len(set(var.ac for var in self.variants if var.ac))
        return  num == 1 or num == 0 
    
    @property
    def all_same_type(self):
        """Return True if all variants have same type, otherwise False."""
        return len(set(var.type for var in self.variants if var.type)) == 1
    
    @property
    def all_same_pos(self):
        """Return True if all variants have same position, otherwise False."""
        return len(set(var.posedit.pos for var in self.variants if var.posedit.pos)) == 1
    
    @property
    def ac(self):
        if self.all_same_ac:
            return self.variants[0].ac
        else:
            return [var.ac for var in self.variants]
    
    @property
    def type(self):
        return self.variants[0].type
    
    @property
    def posedit(self):
        posedits = [var.posedit for var in self.variants]
        return hgvs.posedit.PosEditSet(posedits, self.uncertain)



class CompoundVariant(ComplexVariant):
    """
    represents a compound variant.
    """
    
    separator = ';'
    
    def __init__(self, variants, phases, uncertain=False):
        super(CompoundVariant, self).__init__(variants, phases, uncertain)



class MosaicVariant(ComplexVariant):
    """
    represents a compound variant.
    """
    
    separator = '/'
    
    def __init__(self, variants):
        super(MosaicVariant, self).__init__(variants, None)



class ChimericVariant(ComplexVariant):
    """
    represents a compound variant.
    """
    
    separator = '//'
    
    def __init__(self, variants):
        super(ChimericVariant, self).__init__(variants, None)



## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
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
