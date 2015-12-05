# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import hgvs
import recordtype


class PosEdit(recordtype.recordtype('PosEdit', [('pos', None), ('edit', None), ('uncertain', False)])):
    """
    represents a **simple** variant, consisting of a single position and edit pair
    """
    def __str__(self):
        if not hasattr(self, 'formatting'):
            self.formatting = {}
        
        if self.pos is None:
            rv = str(self.edit.format(self.formatting))
        else:
            rv = '{pos}{edit}'.format(pos=self.pos.format(self.formatting), edit=self.edit.format(self.formatting))
        
        if self.uncertain:
            if self.edit in ['0', '']:
                rv = rv + '?'
            else:
                rv = '(' + rv + ')'
        return rv

    def _set_uncertain(self):
        """sets the uncertain flag to True; used primarily by the HGVS grammar

        :returns: self
        """
        self.uncertain = True
        return self
    
    def format(self, conf=None):
        """Formatting the stringification of PosEdit

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
