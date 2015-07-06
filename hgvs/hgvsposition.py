# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import recordtype


class HGVSPosition(recordtype.recordtype('HGVSPosition', ['ac', 'type', 'pos'])):
    """
    HGVSPosition -- Represent partial HGVS tags that refer to a position without alleles
    
    :param str ac: sequence accession
    :param str type: type of sequence and coordinate
    :param str pos: sequence position
    
    """

    def __str__(self):
        return '{self.ac}:{self.type}.{self.pos}'.format(self=self)

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
