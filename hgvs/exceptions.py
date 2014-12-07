# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

class HGVSError(Exception):
    pass

class HGVSParseError(HGVSError):
    pass

class HGVSInvalidIntervalError(HGVSError):
    pass

class HGVSInvalidVariantError(HGVSError):
    pass

class HGVSValidationError(HGVSError):
    pass

class HGVSDataNotAvailableError(HGVSError):
    pass


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
