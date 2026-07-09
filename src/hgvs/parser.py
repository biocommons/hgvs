# -*- coding: utf-8 -*-
"""Provides parser for HGVS strings and HGVS-related conceptual
components, such as intronic-offset coordiates

"""

import copy
import logging
import re
import warnings

import bioutils.sequences
import ometa.runtime
import parsley
from pyparsing import ParseException

import hgvs.edit

# The following imports are referenced by fully-qualified name in the
# OMeta hgvs grammar.
import hgvs.enums
import hgvs.grammar
import hgvs.hgvsposition
import hgvs.location
import hgvs.posedit
import hgvs.sequencevariant
from hgvs.exceptions import HGVSParseError
from hgvs.generated.hgvs_grammar import createParserClass

#: ``grammar_fn`` value selecting the new pyparsing grammar (``src/hgvs/grammar.py``).
PYPARSING_GRAMMAR = "__pyparsing__"

#: ``grammar_fn`` value selecting the traditional OMeta/parsley grammar.
OMETA_GRAMMAR = "__ometa__"

# Sentinel used to detect whether the caller explicitly passed ``grammar_fn``.
# We can't use ``None`` because we want to distinguish "not specified" (use the
# default backend, no warning) from any explicitly requested value.
_UNSET_GRAMMAR_FN = object()


class _GrammarProxy:
    """Mimics Parsley's grammar(input).rule() calling convention."""

    def __init__(self, grammar_obj, input_string):
        self._grammar_obj = grammar_obj
        self._input = input_string

    def __getattr__(self, rule_name):
        def parse_fn():
            return self._grammar_obj.parse(rule_name, self._input)
        return parse_fn


class _GrammarCallable:
    """Callable that mimics Parsley's grammar API for backwards compatibility.

    Supports: grammar(input).rule(), dir(grammar._grammarClass)
    """

    def __init__(self, grammar_obj):
        self._grammar_obj = grammar_obj
        self._grammarClass = self  # for _expose_rule_functions compatibility

    def __call__(self, input_string):
        return _GrammarProxy(self._grammar_obj, input_string)

    def __dir__(self):
        return ["rule_" + name for name in self._grammar_obj.rules]


class Parser:
    """Provides comprehensive parsing of HGVS variant strings (*i.e.*,
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
    SequenceVariant(ac=NM_01234.5, type=c, posedit=22+1A>T, gene=None)
    >>> v.posedit.pos
    BaseOffsetInterval(start=22+1, end=22+1, uncertain=False)
    >>> i = hp.parse_c_interval("22+1")
    >>> i
    BaseOffsetInterval(start=22+1, end=22+1, uncertain=False)

    The `parse_hgvs_variant` and `parse_c_interval` methods correspond
    to the `hgvs_variant` and `c_interval rules` in the grammar,
    respectively.

    As a convenience, the Parser provides the `parse` method as a
    shorthand for `parse_hgvs_variant`:
    >>> v = hp.parse("NM_01234.5:c.22+1A>T")
    >>> v
    SequenceVariant(ac=NM_01234.5, type=c, posedit=22+1A>T, gene=None)

    Because the methods are generated on-the-fly and depend on the
    grammar that is loaded at runtime, a full list of methods is not
    available in the documentation.  However, the list of
    rules/methods is available via the `rules` instance variable.

    A few notable methods are listed below:

    `parse_hgvs_variant()` parses any valid HGVS string supported by the grammar.

      >>> hp.parse_hgvs_variant("NM_01234.5:c.22+1A>T")
      SequenceVariant(ac=NM_01234.5, type=c, posedit=22+1A>T, gene=None)
      >>> hp.parse_hgvs_variant("NP_012345.6:p.Ala22Trp")
      SequenceVariant(ac=NP_012345.6, type=p, posedit=Ala22Trp, gene=None)

    The `hgvs_variant` rule iteratively attempts parsing using the
    major classes of HGVS variants. For slight improvements in
    efficiency, those rules may be invoked directly:

      >>> hp.parse_p_variant("NP_012345.6:p.Ala22Trp")
      SequenceVariant(ac=NP_012345.6, type=p, posedit=Ala22Trp, gene=None)

    Similarly, components of the underlying structure may be parsed
    directly as well:

      >>> hp.parse_c_posedit("22+1A>T")
      PosEdit(pos=22+1, edit=A>T, uncertain=False)
      >>> hp.parse_c_interval("22+1")
      BaseOffsetInterval(start=22+1, end=22+1, uncertain=False)

    Two grammar implementations are available, selected via ``grammar_fn``:

    * By default, the traditional OMeta/parsley grammar is used, preserving
      historical behavior.
    * ``Parser(grammar_fn="__pyparsing__")`` selects the newer pyparsing
      grammar (see ``hgvs.grammar``). This is expected to become the default
      in a future major release, at which point OMeta/parsley will be
      deprecated and eventually removed.
    * ``Parser(grammar_fn="__ometa__")`` explicitly selects the OMeta/parsley
      grammar (equivalent to the current default).

    Passing any other value (a path to a custom OMeta grammar file) is
    deprecated and will be removed in a future version.

    """

    def __init__(self, grammar_fn=_UNSET_GRAMMAR_FN, expose_all_rules=False):
        self._logger = logging.getLogger(__name__)
        if grammar_fn is _UNSET_GRAMMAR_FN or grammar_fn == OMETA_GRAMMAR:
            # Default (1.x): the traditional OMeta/parsley grammar.
            self._grammar = self._make_ometa_grammar()
        elif grammar_fn == PYPARSING_GRAMMAR:
            self._grammar_obj = hgvs.grammar.HGVSGrammar()
            self._grammar = _GrammarCallable(self._grammar_obj)
        else:
            # Deprecated escape hatch: a path to a custom OMeta grammar file.
            warnings.warn(
                "Passing a grammar filename to Parser(grammar_fn=...) is deprecated "
                "and will be removed in a future version. Use grammar_fn='__pyparsing__' "
                "to select the pyparsing grammar or grammar_fn='__ometa__' for the "
                "traditional grammar.",
                DeprecationWarning,
                stacklevel=2,
            )
            self._grammar = self._make_ometa_grammar(grammar_fn)
        self._expose_rule_functions(expose_all_rules)

    @staticmethod
    def _make_ometa_grammar(grammar_fn=None):
        """Build the traditional OMeta/parsley grammar.

        With ``grammar_fn=None`` the bundled, pre-generated grammar is used.
        Otherwise ``grammar_fn`` is read as a custom OMeta grammar file.
        """
        bindings = {"hgvs": hgvs, "bioutils": bioutils, "copy": copy}
        if grammar_fn is None:
            return parsley.wrapGrammar(
                createParserClass(ometa.runtime.OMetaGrammarBase, bindings)
            )
        with open(grammar_fn, "r") as grammar_file:
            return parsley.makeGrammar(grammar_file.read(), bindings)

    def parse(self, v) -> hgvs.sequencevariant.SequenceVariant:
        """parse HGVS variant `v`, returning a SequenceVariant

        :param str v: an HGVS-formatted variant as a string
        :rtype: SequenceVariant

        """
        return self.parse_hgvs_variant(v)

    def _expose_rule_functions(self, expose_all_rules=False):
        """add parse functions for public grammar rules

        Defines a function for each public grammar rule, based on
        introspecting the grammar. For example, the `c_interval` rule
        is exposed as a method `parse_c_interval` and used like this::

          Parser.parse_c_interval('26+2_57-3') -> Interval(...)

        """

        def make_parse_rule_function(rule_name):
            "builds a wrapper function that parses a string with the specified rule"

            def rule_fxn(s):
                try:
                    return self._grammar(s).__getattr__(rule_name)()
                except ometa.runtime.ParseError as exc:
                    # OMeta/parsley backend
                    raise HGVSParseError(
                        "{s}: char {exc.position}: {reason}".format(
                            s=s, exc=exc, reason=exc.formatReason()
                        )
                    )
                except ParseException as exc:
                    # pyparsing backend
                    raise HGVSParseError(
                        "{s}: char {pos}: {msg}".format(
                            s=s, pos=exc.loc, msg=exc.msg
                        )
                    )

            rule_fxn.__doc__ = "parse string s using `%s' rule" % rule_name
            return rule_fxn

        exposed_rule_re = re.compile(
            r"hgvs_(variant|position)|(c|g|m|n|p|r)"
            r"_(edit|hgvs_position|interval|pos|posedit|variant)"
        )
        exposed_rules = [
            m.replace("rule_", "")
            for m in dir(self._grammar._grammarClass)
            if m.startswith("rule_")
        ]
        if not expose_all_rules:
            exposed_rules = [
                rule_name for rule_name in exposed_rules if exposed_rule_re.match(rule_name)
            ]
        for rule_name in exposed_rules:
            att_name = "parse_" + rule_name
            rule_fxn = make_parse_rule_function(rule_name)
            self.__setattr__(att_name, rule_fxn)
        self._logger.debug(
            "Exposed {n} rules ({rules})".format(
                n=len(exposed_rules), rules=", ".join(exposed_rules)
            )
        )


# <LICENSE>
# Copyright 2018 HGVS Contributors (https://github.com/biocommons/hgvs)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
