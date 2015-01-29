# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import re


def build_tx_cigar(exons, strand):
    cigarelem_re = re.compile('\d+[=DIMNX]')
    def _reverse_cigar(c):
        return ''.join(reversed(cigarelem_re.findall(c)))

    if len(exons) == 0:
        return None

    tx_cigar = [exons[0]['cigar']]    # exon 1
    for i in range(1, len(exons)):    # and intron + exon pairs thereafter
        cigar = exons[i]['cigar']
        if strand == -1:
            cigar = _reverse_cigar(cigar)
        tx_cigar += [str(exons[i]['alt_start_i'] - exons[i-1]['alt_end_i']) + 'N',
                     cigar]
    
    tx_cigar_str = ''.join(tx_cigar)
    
    return tx_cigar_str





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
