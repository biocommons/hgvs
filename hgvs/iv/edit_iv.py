__doc__ = """
hgvs.edit_iv -- representation of edit operations in Invitae-specific HGVS variants

These are non-standard forms used internally

"""

import recordtype

from hgvs.edit import Edit

class NACopy(Edit, recordtype.recordtype('NACopy', ['copy', ('uncertain', False)])):

    def __str__(self):
        s = 'copy{}'.format(self.copy)
        return '('+s+')' if self.uncertain else s

    def set_uncertain(self):
        self.uncertain = True
        return self


class NADupN(Edit, recordtype.recordtype('NADupN', ['n', ('uncertain', False)])):

    def __str__(self):
        s = 'dup{}'.format(self.n)
        return '('+s+')' if self.uncertain else s

    def set_uncertain(self):
        self.uncertain = True
        return self

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
