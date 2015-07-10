# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

from pkg_resources import resource_filename

import bioutils.sequences
import ometa.runtime
import parsley

from .exceptions import HGVSParseError
import hgvs.edit
import hgvs.hgvsposition
import hgvs.location
import hgvs.posedit
import hgvs.variant


class Parser(object):
    """Provides comprehensive parsing of HGVS varaint strings (*i.e.*,
    variants represented according to the Human Genome Variation
    Society recommendations) into Python representations.  The class
    wraps a Parsing Expression Grammar, exposing rules of that grammar
    as methods (prefixed with `parse_`) that parse an input string
    according to the rule.  The class exposes all rules, so that it's
    possible to parse both full variant representations as well as
    components, like so:
    
    >>> hp = Parser()
    >>> v = hp.parse_hgvs_variant("NM_01234.5:c.22+1A>T")
    >>> v
    SequenceVariant(ac=NM_01234.5, type=c, posedit=22+1A>T)
    >>> v.posedit.pos
    Interval(start=22+1, end=22+1, uncertain=False)
    >>> i = hp.parse_c_interval('22+1')
    >>> i
    Interval(start=22+1, end=22+1, uncertain=False)

    The `parse_hgvs_variant` and `parse_c_interval` methods correspond
    to the `hgvs_variant` and `c_interval rules` in the grammar,
    respectively.

    Because the methods are generated on-the-fly and depend on the
    grammar that is loaded at runtime, a full list of methods is not
    available in the documentation.  However, the list of
    rules/methods is available via the `rules` instance variable.

    A few notable methods are listed below:

    `parse_hgvs_variant()` parses any valid HGVS string supported by the grammar.

      >>> hp.parse_hgvs_variant("NM_01234.5:c.22+1A>T")
      SequenceVariant(ac=NM_01234.5, type=c, posedit=22+1A>T)
      >>> hp.parse_hgvs_variant('NP_012345.6:p.Ala22Trp')
      SequenceVariant(ac=NP_012345.6, type=p, posedit=Ala22Trp)

    The `hgvs_variant` rule iteratively attempts parsing using the
    major classes of HGVS variants. For slight improvements in
    efficiency, those rules may be invoked directly:

      >>> hp.parse_p_variant('NP_012345.6:p.Ala22Trp')
      SequenceVariant(ac=NP_012345.6, type=p, posedit=Ala22Trp)

    Similarly, components of the underlying structure may be parsed
    directly as well:

      >>> hp.parse_c_posedit('22+1A>T')
      PosEdit(pos=22+1, edit=A>T, uncertain=False)
      >>> hp.parse_c_interval('22+1')
      Interval(start=22+1, end=22+1, uncertain=False)

    """

    __default_grammar_fn = resource_filename(__name__, '_data/hgvs.pymeta')

    def __init__(self, grammar_fn=__default_grammar_fn):
        self._grammar_fn = grammar_fn
        self._grammar = parsley.makeGrammar(open(grammar_fn, 'r').read(), {'hgvs': hgvs, 'bioutils': bioutils})

        # define function attributes for each grammar rule, prefixed with 'parse_'
        # e.g., Parser.parse_c_interval('26+2_57-3') -> Interval(...)
        #TODO: exclude built-in rules
        self.rules = [m.replace('rule_', '') for m in dir(self._grammar._grammarClass) if m.startswith('rule_')]
        for rule_name in self.rules:
            att_name = 'parse_' + rule_name
            rule_fxn = self.__make_parse_rule_function(rule_name)
            self.__setattr__(att_name, rule_fxn)

    def __make_parse_rule_function(self, rule_name):
        # http://docs.python.org/2/reference/datamodel.html#object.__getattr__
        """
        This function returns a function that takes a string and returns the parsing result.
        """

        def rule_fxn(s):
            try:
                return self._grammar(s).__getattr__(rule_name)()
            except ometa.runtime.ParseError as exc:
                raise HGVSParseError(
                    "{s}: char {exc.position}: {reason}".format(s=s,
                                                                exc=exc,
                                                                reason=exc.formatReason()))

        rule_fxn.func_doc = "parse string s using `%s' rule" % rule_name
        return rule_fxn

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
