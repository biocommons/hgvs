# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
"""
hgvs.formatter
"""

import copy
import logging

import hgvs
import hgvs.parser
import hgvs.variantmapper
import hgvs.dataproviders.uta

_logger = logging.getLogger(__name__)



class Formatter(object):
    """Perform variant formatting
    """

    def __init__(self, hdp=None,
                 fill_ref=hgvs.global_config.formatting.fill_ref,
                 ):
        """Initialize and configure the formatter

        :param hdp: HGVS Data Provider Interface-compliant instance (see :class:`hgvs.dataproviders.interface.Interface`)
        :param fill_ref: fill in reference allele for eidts
        """
        self.hdp = hdp
        self.fill_ref = fill_ref
        if self.fill_ref:
            self.hm = hgvs.variantmapper.VariantMapper(self.hdp)
    
    
    def format(self, var):
        var_i = copy.deepcopy(var)
        self._format_ref(var_i)
        return var_i
    
    def _format_ref(self, var):
        if self.fill_ref:
            if var.posedit.edit.type in ['del', 'delins', 'identity', 'dup', 'inv'] and var.posedit.edit.ref_s is None:
                self.hm._replace_reference_for_sequence_variant(var)
        else:
            if var.posedit.edit.type in ['del', 'delins', 'identity', 'dup', 'inv']:
                var.posedit.edit.ref = ''
        return var



if __name__ == '__main__':
    hgvsparser = hgvs.parser.Parser()
    var = hgvsparser.parse_hgvs_variant('NM_001166478.1:c.61del')
    hdp = hgvs.dataproviders.uta.connect()
    hf  = Formatter(hdp)
    res = hf.format(var)
    print(str(var) + '    =>    ' + str(res))


## <LICENSE>
## Copyright 2015 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
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
