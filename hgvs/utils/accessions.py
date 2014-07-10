import warnings

NC_to_chr_dict = {
    'NC_000001.10':  '1', 'NC_000002.11':  '2', 'NC_000003.11': '3',
    'NC_000004.11':  '4', 'NC_000005.9' :  '5', 'NC_000006.11': '6',
    'NC_000007.13':  '7', 'NC_000008.10':  '8', 'NC_000009.11': '9',
    'NC_000010.10': '10', 'NC_000011.9' : '11', 'NC_000012.11': '12',
    'NC_000013.10': '13', 'NC_000014.8' : '14', 'NC_000015.9' : '15',
    'NC_000016.9' : '16', 'NC_000017.10': '17', 'NC_000018.9' : '18',
    'NC_000019.9' : '19', 'NC_000020.10': '20', 'NC_000021.8' : '21',
    'NC_000022.10': '22', 'NC_000023.10':  'X', 'NC_000024.9' : 'Y',
    }

chr_to_NC_dict = dict([ (v,k) for k,v in NC_to_chr_dict.iteritems() ])

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
