import warnings

import recordtype

class SequenceVariant( recordtype.recordtype('Variant', ['ac','type','posedit']) ):
    """
    represents a basic HGVS variant.  The only requirement is that each
    component can be stringified; for example, passing pos as either a string
    or an hgvs.location.CDSInterval (for example) are both intended uses
    """
    
    @property
    def seqref(self):
        warnings.warn('seqref is deprecated; use ac instead',DeprecationWarning,stacklevel=2)
        return self.ac

    def __str__(self):
        return '{self.ac}:{self.type}.{self.posedit}'.format(self=self)

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
