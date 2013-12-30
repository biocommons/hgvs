import recordtype

from hgvs.exceptions import HGVSError

class PosEdit( recordtype.recordtype( 'PosEdit', [('pos',None),('edit',None),('uncertain',False)] )):
    """
    represents a **simple** variant
    """
    
    def __str__(self):
        rv = str(self.edit) if self.pos is None else '{self.pos}{self.edit}'.format(self=self)
        if self.uncertain:
            if self.edit in ['0','']:
                rv = rv + '?'
            else:
                rv = '(' + rv + ')'
        return rv

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
