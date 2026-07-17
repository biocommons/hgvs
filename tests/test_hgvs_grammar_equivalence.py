# -*- coding: utf-8 -*-
"""Tests that the pyparsing and OMeta/parsley backends behave identically.

The pyparsing grammar is intended to replace OMeta/parsley (issue #505). While
both are shipped, the OMeta backend is the reference implementation: any
difference in what is accepted, rejected, or returned is a bug in the new
grammar unless it is a deliberate, documented change.

The one intentional difference is leading whitespace -- see
TestWhitespaceHandling below.
"""

import os
import unittest

import pyparsing as pp

import hgvs.parser
from hgvs.exceptions import HGVSParseError

BACKENDS = (hgvs.parser.OMETA_GRAMMAR, hgvs.parser.PYPARSING_GRAMMAR)

_GAUNTLET_FN = os.path.join(os.path.dirname(__file__), "data", "gauntlet")


def _read_gauntlet():
    with open(_GAUNTLET_FN, "r") as f:
        return [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]


def _parse_or_error(parser, variant):
    """Return the stringified parse, or None if the variant was rejected."""
    try:
        return str(parser.parse_hgvs_variant(variant))
    except HGVSParseError:
        return None


class TestBackendEquivalence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ometa = hgvs.parser.Parser(grammar_fn=hgvs.parser.OMETA_GRAMMAR)
        cls.pyparsing = hgvs.parser.Parser(grammar_fn=hgvs.parser.PYPARSING_GRAMMAR)

    def test_backends_agree_on_gauntlet(self):
        """both backends must return identical results for every gauntlet variant"""

        variants = _read_gauntlet()
        self.assertGreater(len(variants), 0, "gauntlet is empty; test would be vacuous")

        for variant in variants:
            with self.subTest(variant=variant):
                self.assertEqual(
                    _parse_or_error(self.ometa, variant),
                    _parse_or_error(self.pyparsing, variant),
                    "backends disagree on {!r}".format(variant),
                )

    def test_exposed_rules_match(self):
        """the public parse_* API must be identical between backends"""

        def parse_methods(parser):
            return set(a for a in dir(parser) if a.startswith("parse_"))

        self.assertEqual(parse_methods(self.ometa), parse_methods(self.pyparsing))


class TestNoGlobalPyparsingState(unittest.TestCase):
    """pyparsing's whitespace/packrat settings are process-global.

    hgvs must not leak them: doing so silently changes the behavior of any
    other pyparsing grammar in the same process. cf. hgvs.grammar.
    """

    def test_building_grammar_does_not_leak_whitespace_default(self):
        before = pp.ParserElement.DEFAULT_WHITE_CHARS
        hgvs.parser.Parser(grammar_fn=hgvs.parser.PYPARSING_GRAMMAR)
        self.assertEqual(pp.ParserElement.DEFAULT_WHITE_CHARS, before)

    def test_third_party_grammar_unaffected(self):
        hgvs.parser.Parser(grammar_fn=hgvs.parser.PYPARSING_GRAMMAR)
        expr = pp.Word(pp.alphas) + pp.Word(pp.nums)
        self.assertEqual(expr.parse_string("abc 123").as_list(), ["abc", "123"])

    def test_grammar_is_strict_regardless_of_ambient_whitespace(self):
        """a hostile global whitespace default must not loosen the grammar"""

        prev = pp.ParserElement.DEFAULT_WHITE_CHARS
        try:
            pp.ParserElement.set_default_whitespace_chars(" \t\n")
            parser = hgvs.parser.Parser(grammar_fn=hgvs.parser.PYPARSING_GRAMMAR)
            for variant in (
                "NM_01234.5 :c.22+1A>T",
                "NM_01234.5: c.22+1A>T",
                "NM_01234.5:c.22 +1A>T",
            ):
                with self.subTest(variant=variant):
                    with self.assertRaises(HGVSParseError):
                        parser.parse_hgvs_variant(variant)
        finally:
            pp.ParserElement.set_default_whitespace_chars(prev)


class TestWhitespaceHandling(unittest.TestCase):
    """Whitespace is the one intentional behavior change (issue #505).

    OMeta silently skips leading spaces/tabs; the pyparsing grammar rejects
    them. HGVS strings contain no whitespace, so the stricter reading is
    correct. Both backends reject trailing whitespace and leading newlines.
    """

    VARIANT = "NM_01234.5:c.22+1A>T"

    def test_ometa_accepts_leading_space_pyparsing_rejects(self):
        ometa = hgvs.parser.Parser(grammar_fn=hgvs.parser.OMETA_GRAMMAR)
        pyparsing = hgvs.parser.Parser(grammar_fn=hgvs.parser.PYPARSING_GRAMMAR)
        for padded in (" " + self.VARIANT, "\t" + self.VARIANT, "  " + self.VARIANT):
            with self.subTest(variant=padded):
                self.assertEqual(str(ometa.parse_hgvs_variant(padded)), self.VARIANT)
                with self.assertRaises(HGVSParseError):
                    pyparsing.parse_hgvs_variant(padded)

    def test_both_backends_reject_trailing_whitespace(self):
        for backend in BACKENDS:
            parser = hgvs.parser.Parser(grammar_fn=backend)
            for padded in (self.VARIANT + " ", self.VARIANT + "\t", "\n" + self.VARIANT):
                with self.subTest(backend=backend, variant=padded):
                    with self.assertRaises(HGVSParseError):
                        parser.parse_hgvs_variant(padded)


if __name__ == "__main__":
    unittest.main()

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
