"""pyparsing backend for the HGVS parser (see hgvs.parsers.pyparsing_grammar)."""

from pyparsing import ParseBaseException

import hgvs.parsers.pyparsing_grammar
from hgvs.exceptions import HGVSParseError
from hgvs.parsers.base import Parser


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
        # Name mirrors Parsley's API, which _expose_rule_functions introspects.
        self._grammarClass = self  # NOSONAR: S116 - matches the Parsley attribute name

    def __call__(self, input_string):
        return _GrammarProxy(self._grammar_obj, input_string)

    def __dir__(self):
        return ["rule_" + name for name in self._grammar_obj.rules]


class PyParsingParser(Parser):
    """The newer pyparsing grammar backend. See ``Parser`` for the public API."""

    # ParseBaseException rather than ParseException: ParseFatalException and
    # ParseSyntaxException are siblings of ParseException, not subclasses of it.
    _parse_error_cls = ParseBaseException

    def _build_grammar(self, grammar_fn):
        self._grammar_obj = hgvs.parsers.pyparsing_grammar.HGVSGrammar()
        return _GrammarCallable(self._grammar_obj)

    def _translate_parse_error(self, s, exc):
        return HGVSParseError(f"{s}: char {exc.loc}: {exc.msg}")


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
