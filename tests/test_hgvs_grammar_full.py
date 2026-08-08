import csv
import os
import pprint
import re
import unittest
from importlib.resources import files as resources_files

import hgvs.parsers

#
# Tests of the grammar
#
# Code takes a tab-delimited text file of the form:
#
#   Func    Test    Valid   InType  Expected
#   pm      -+      True    string
#   pm      *       False   one
#   num     1|+1    True    list    1|1
#
# Headers are defined as follows:
# Func: function name to call in the grammar
# Test: item(s) to test
# Valid: if the input is expected to be valid (True or False)
# InType: 3 type:
# - one: input is a single value
# - string: input is a string; test each character in the string separately
# - list: input is a list delimited by a pipe character ("|")
# Expected: expected result (if stringifying input does not return the same answer, e,g. "+1" -> "1")
# - if expected is left blank, then it is assumed that stringifying the parsed input returns the same answer.
#


#: Both grammar backends are exercised: the OMeta/parsley grammar is the 1.x
#: default, and pyparsing is slated to become the default in 2.x.
BACKENDS = (hgvs.parsers.OMETA_GRAMMAR, hgvs.parsers.PYPARSING_GRAMMAR)


class TestGrammarFull(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p = hgvs.parsers.Parser()
        cls.grammar = cls.p._grammar
        cls._test_fn = os.path.join(os.path.dirname(__file__), "data", "grammar_test.tsv")

    def test_parser_test_completeness(self):
        """ensure that all rules in grammar have tests"""

        from hgvs.parsers.pyparsing_grammar import HGVSGrammar

        grammar_rules = set(HGVSGrammar().rules.keys())

        with open(self._test_fn) as f:
            reader = csv.DictReader(f, delimiter="\t")
            test_rules = set(row["Func"] for row in reader)

        untested_rules = grammar_rules - test_rules

        self.assertTrue(len(untested_rules) == 0, f"untested rules: {untested_rules}")

    def test_grammar_rules_match_between_backends(self):
        """both backends must expose exactly the same rule names"""

        from hgvs.parsers.pyparsing_grammar import HGVSGrammar

        grammar_rule_re = re.compile(r"^(\w+)")
        grammar_fn = resources_files("hgvs.parsers") / "_data" / "hgvs.pymeta"
        with open(grammar_fn) as f:
            ometa_rules = set(r.group(1) for r in filter(None, map(grammar_rule_re.match, f)))
        pyparsing_rules = set(HGVSGrammar().rules.keys())

        self.assertEqual(
            ometa_rules,
            pyparsing_rules,
            f"rules differ between backends: only in OMeta: {sorted(ometa_rules - pyparsing_rules)}; only in pyparsing: {sorted(pyparsing_rules - ometa_rules)}",
        )

    def test_parser_grammar(self):
        for backend in BACKENDS:
            with self.subTest(backend=backend):
                self._check_grammar(hgvs.parsers.Parser(grammar_fn=backend))

    def _check_grammar(self, parser):
        with open(self._test_fn) as f:
            reader = csv.DictReader(f, delimiter="\t")

            fail_cases = []

            for row in reader:
                if row["Func"].startswith("#"):
                    continue

                # setup input
                inputs = self._split_inputs(row["Test"], row["InType"])
                expected_results = (
                    self._split_inputs(row["Expected"], row["InType"])
                    if row["Expected"]
                    else inputs
                )
                expected_map = dict(zip(inputs, expected_results, strict=False))
                # step through each item and check
                is_valid = True if row["Valid"].lower() == "true" else False

                for key in expected_map:
                    expected_result = str(expected_map[key]).replace("u'", "'")
                    function_to_test = getattr(parser._grammar(key), row["Func"])
                    row_str = "{}\t{}\t{}\t{}\t{}".format(
                        row["Func"], key, row["Valid"], "one", expected_result
                    )
                    try:
                        actual_result = str(function_to_test()).replace("u'", "'")
                        if not is_valid or (expected_result != actual_result):
                            print(f"expected: {expected_result} actual:{actual_result}")
                            fail_cases.append(row_str)
                    except Exception as e:
                        if is_valid:
                            print(f"expected: {expected_result} Exception: {e}")
                            fail_cases.append(row_str)

        # everything should have passed - report whatever failed
        self.assertTrue(len(fail_cases) == 0, pprint.pprint(fail_cases))

    def _split_inputs(self, in_string, intype):
        DELIM = "|"
        if intype == "list":
            inputs = in_string.split(DELIM)
        elif intype == "string":
            inputs = list(in_string)
        elif intype == "one":
            inputs = [in_string]
        else:
            assert False, f"shouldn't be here (intype = {intype})"
        inputs = [x if x != "None" else None for x in inputs]
        return inputs


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
