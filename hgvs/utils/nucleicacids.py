# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import string

complement_transtable_str = string.maketrans('ACGT','TGCA')
complement_transtable_uni = dict(zip(map(ord,u'ACGT'),u'TGCA'))
def reverse_complement(s):
    def _rc_str(s):
        return ''.join(reversed(s.translate(complement_transtable_str)))
    def _rc_uni(s):
        return ''.join(reversed(s.translate(complement_transtable_uni)))
    return None if s is None else _rc_uni(s) if isinstance(s,unicode) else _rc_str(s)




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
