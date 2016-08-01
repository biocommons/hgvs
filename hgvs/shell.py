# -*- coding: utf-8 -*-

"""start IPython shell with hgvs initialized. Intended to be used for
experimenting, debugging, and generating bug reports."""

from __future__ import absolute_import, division, print_function, unicode_literals

import logging

import IPython

import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.parser
import hgvs.variantmapper
import hgvs.validator

logging.basicConfig(level=logging.INFO)

alt_aln_method = 'splign'
hgvs_g = 'NC_000007.13:g.36561662C>T'
hgvs_c = 'NM_001637.3:c.1582G>A'
hgvs_p = 'NP_001628.1:p.(Gly528Arg)'


def shell():
    hp = hgvsparser = hgvs.parser.Parser()

    hdp = hgvs.dataproviders.uta.connect()
    vm = variantmapper = hgvs.variantmapper.VariantMapper(hdp)
    evm = easyvariantmapper = hgvs.variantmapper.EasyVariantMapper(
        hdp, primary_assembly='GRCh37', alt_aln_method=alt_aln_method)
    hv = hgvs.validator.Validator(hdp)
    hn = hgvs.normalizer.Normalizer(hdp)

    var_g = hgvsparser.parse_hgvs_variant(hgvs_g)
    var_c = hgvsparser.parse_hgvs_variant(hgvs_c)
    var_p = hgvsparser.parse_hgvs_variant(hgvs_p)

    IPython.embed()



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
