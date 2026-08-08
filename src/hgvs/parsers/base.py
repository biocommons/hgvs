"""Base class for HGVS parser backends.

``Parser`` is the public entry point. It never constructs itself directly;
``__new__`` dispatches to a backend-specific subclass (``ParsleyParser`` or
``PyParsingParser``) based on ``grammar_fn``, so that each backend's
dependencies and parsing internals stay isolated in their own module and are
only imported when that backend is actually selected.
"""

import logging
import re

import hgvs.sequencevariant

#: ``grammar_fn`` value selecting the pyparsing grammar (``PyParsingParser``).
PYPARSING_GRAMMAR = "__pyparsing__"

#: ``grammar_fn`` value selecting the traditional OMeta/parsley grammar (``ParsleyParser``).
OMETA_GRAMMAR = "__ometa__"

# Sentinel used to detect whether the caller explicitly passed ``grammar_fn``.
# We can't use ``None`` because we want to distinguish "not specified" (use the
# default backend, no warning) from any explicitly requested value.
_UNSET_GRAMMAR_FN = object()

_EXPOSED_RULE_RE = re.compile(
    r"hgvs_(variant|position)|(c|g|m|n|p|r)" r"_(edit|hgvs_position|interval|pos|posedit|variant)"
)


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
      historical behavior. ``Parser()`` returns a ``ParsleyParser`` instance.
    * ``Parser(grammar_fn="__pyparsing__")`` selects the newer pyparsing
      grammar, returning a ``PyParsingParser`` instance. This is expected to
      become the default in a future major release, at which point
      OMeta/parsley will be deprecated and eventually removed.
    * ``Parser(grammar_fn="__ometa__")`` explicitly selects the OMeta/parsley
      grammar (equivalent to the current default).

    Passing any other value (a path to a custom OMeta grammar file) is
    deprecated and will be removed in a future version.

    The two backends are intended to be interchangeable, and are tested for
    equivalence (see tests/test_hgvs_grammar_equivalence.py). There is one
    deliberate difference: the OMeta grammar silently skips leading spaces and
    tabs, so ``Parser().parse_hgvs_variant(" NM_01234.5:c.22+1A>T")`` succeeds,
    whereas the pyparsing grammar rejects it. HGVS strings contain no
    whitespace, so the stricter reading is intended. Both backends reject
    trailing whitespace and leading newlines. Callers relying on the older,
    laxer behavior should strip the input themselves.

    """

    def __new__(cls, grammar_fn=_UNSET_GRAMMAR_FN, **kwargs):
        if cls is not Parser:
            # A backend subclass (or a further subclass of one) was
            # constructed directly -- don't re-dispatch.
            return super().__new__(cls)

        # Local imports: each backend module is only imported when that
        # backend is actually selected, so callers who use only one backend
        # never pay the import cost of the other (parsley+ometa, or pyparsing).
        if grammar_fn is _UNSET_GRAMMAR_FN or grammar_fn == OMETA_GRAMMAR:
            from hgvs.parsers.parsley_parser import ParsleyParser

            return super().__new__(ParsleyParser)
        if grammar_fn == PYPARSING_GRAMMAR:
            from hgvs.parsers.pyparsing_parser import PyParsingParser

            return super().__new__(PyParsingParser)
        # Deprecated escape hatch: a path to a custom OMeta grammar file.
        from hgvs.parsers.parsley_parser import ParsleyParser

        return super().__new__(ParsleyParser)

    def __init__(self, grammar_fn=_UNSET_GRAMMAR_FN, expose_all_rules=False):
        self._logger = logging.getLogger(__name__)
        self._grammar = self._build_grammar(grammar_fn)
        self._expose_rule_functions(expose_all_rules)

    def _build_grammar(self, grammar_fn):
        """Build and return the backend-specific grammar object. Implemented by subclasses."""
        raise NotImplementedError

    #: The backend-native exception type raised on a parse failure. Set by subclasses.
    _parse_error_cls = None

    def _translate_parse_error(self, s, exc):
        """Translate a backend-native parse exception into HGVSParseError. Implemented by subclasses."""
        raise NotImplementedError

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
                except self._parse_error_cls as exc:
                    raise self._translate_parse_error(s, exc) from exc

            rule_fxn.__doc__ = f"parse string s using `{rule_name}' rule"
            return rule_fxn

        exposed_rules = [
            m.replace("rule_", "")
            for m in dir(self._grammar._grammarClass)
            if m.startswith("rule_")
        ]
        if not expose_all_rules:
            exposed_rules = [
                rule_name for rule_name in exposed_rules if _EXPOSED_RULE_RE.match(rule_name)
            ]
        for rule_name in exposed_rules:
            att_name = "parse_" + rule_name
            rule_fxn = make_parse_rule_function(rule_name)
            self.__setattr__(att_name, rule_fxn)
        self._logger.debug("Exposed %d rules (%s)", len(exposed_rules), ", ".join(exposed_rules))


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
