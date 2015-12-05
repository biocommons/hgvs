# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

from .decorators.deprecated import deprecated

import recordtype

import hgvs.variantmapper


class SequenceVariant(recordtype.recordtype('SequenceVariant', ['ac', 'type', 'posedit'])):
    """
    represents a basic HGVS variant.  The only requirement is that each
    component can be stringified; for example, passing pos as either a string
    or an hgvs.location.CDSInterval (for example) are both intended uses
    """

    def __str__(self):
        if not hasattr(self, 'formatting'):
            self.formatting = {}
        
        if self.ac is not None:
            return '{ac}:{type}.{posedit}'.format(ac=self.ac, type=self.type, posedit=self.posedit.format(self.formatting))
        else:
            return '{type}.{posedit}'.format(type=self.type, posedit=self.posedit.format(self.formatting))

    
    def format(self, conf=None):
        """Formatting the stringification of sequence variants

        :param conf: a dict comprises formatting options. None is to use global settings.
        formatting configuration options:
            p_3_letter: use 1-letter or 3-letter amino acid representation for p. variants.
            p_term_asterisk: use * or Ter to represent stop-codon gain for p. variants.
        """
        self.formatting = {}
        if conf:
            for option in conf:
                self.formatting[option] = conf[option]
        return str(self)


    def fill_ref(self, hdp):
        hm = hgvs.variantmapper.VariantMapper(hdp)
        if self.posedit.edit.ref_s is None:
            hm._replace_reference(self)
        return self


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
