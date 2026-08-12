"""OMeta/parsley backend for the HGVS parser (the historical, 1.x default)."""

import copy
import warnings
from pathlib import Path

import bioutils
import ometa.runtime
import parsley

import hgvs
import hgvs.edit

# The following imports are referenced by fully-qualified name in the
# OMeta hgvs grammar.
import hgvs.enums
import hgvs.hgvsposition
import hgvs.location
import hgvs.posedit
import hgvs.sequencevariant
from hgvs.exceptions import HGVSParseError
from hgvs.parsers.base import _UNSET_GRAMMAR_FN, OMETA_GRAMMAR, Parser
from hgvs.parsers.generated.hgvs_grammar import createParserClass


class ParsleyParser(Parser):
    """The traditional OMeta/parsley grammar backend. See ``Parser`` for the public API."""

    _parse_error_cls = ometa.runtime.ParseError

    def _build_grammar(self, grammar_fn):
        if grammar_fn is _UNSET_GRAMMAR_FN or grammar_fn == OMETA_GRAMMAR:
            # Default (1.x): the bundled, pre-generated grammar.
            bindings = {"hgvs": hgvs, "bioutils": bioutils, "copy": copy}
            return parsley.wrapGrammar(createParserClass(ometa.runtime.OMetaGrammarBase, bindings))
        # Deprecated escape hatch: a path to a custom OMeta grammar file.
        warnings.warn(
            "Passing a grammar filename to Parser(grammar_fn=...) is deprecated "
            "and will be removed in a future version. Use grammar_fn='__pyparsing__' "
            "to select the pyparsing grammar or grammar_fn='__ometa__' for the "
            "traditional grammar.",
            DeprecationWarning,
            stacklevel=3,
        )
        bindings = {"hgvs": hgvs, "bioutils": bioutils, "copy": copy}
        with Path(grammar_fn).open() as grammar_file:
            return parsley.makeGrammar(grammar_file.read(), bindings)

    def _translate_parse_error(self, s, exc):
        return HGVSParseError(f"{s}: char {exc.position}: {exc.formatReason()}")


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
