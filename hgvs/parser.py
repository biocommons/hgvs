import os
from pkg_resources import resource_filename

import parsley

import hgvs.edit
import hgvs.hgvsposition
import hgvs.location
import hgvs.posedit
import hgvs.variant
import hgvs.utils

from hgvs.exceptions import *

class Parser(object):
    __default_grammar_fn = resource_filename(__name__, 'data/hgvs.pymeta')

    def __init__(self,grammar_fn=__default_grammar_fn):
        self._grammar_fn = grammar_fn
        self._grammar = parsley.makeGrammar(open(grammar_fn,'r').read(),{'hgvs': hgvs})

        # define function attributes for each grammar rule, prefixed with 'parse_'
        # e.g., Parser.parse_c_interval('26+2_57-3') -> Interval(...)
        rules = [ m.replace('rule_','')
                  for m in dir(self._grammar._grammarClass)
                  if m.startswith('rule_') ]
        for rule_name in rules:
            att_name = 'parse_' + rule_name
            rule_fxn = self.__make_parse_rule_function(rule_name)
            self.__setattr__(att_name,rule_fxn)

    def __make_parse_rule_function(self,rule_name):
        # http://docs.python.org/2/reference/datamodel.html#object.__getattr__
        """
        This function returns a function that takes a string and returns the parsing result.
        """
        rule_fxn = lambda s: self._grammar(s).__getattr__(rule_name)()
        rule_fxn.func_doc = "parse string s using `%s' rule" % rule_name
        return rule_fxn

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
