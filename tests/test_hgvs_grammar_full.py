# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

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

import os
import pkg_resources
import pprint
import re
import unittest
from sys import version_info

if version_info < (3, ):
    import unicodecsv as csv
else:
    import csv

import hgvs.parser
from six.moves import map
import six


class TestGrammarFull(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p = hgvs.parser.Parser()
        cls.grammar = cls.p._grammar
        cls._test_fn = os.path.join(os.path.dirname(__file__), "data", "grammar_test.tsv")

    def test_parser_test_completeness(self):
        """ensure that all rules in grammar have tests"""

        grammar_rule_re = re.compile(r"^(\w+)")
        grammar_fn = pkg_resources.resource_filename("hgvs", "_data/hgvs.pymeta")
        with open(grammar_fn, "r") as f:
            grammar_rules = set(r.group(1) for r in filter(None, map(grammar_rule_re.match, f)))

        with open(self._test_fn, "r") as f:
            reader = csv.DictReader(f, delimiter=str("\t"))
            test_rules = set(row["Func"] for row in reader)

        untested_rules = grammar_rules - test_rules

        self.assertTrue(len(untested_rules) == 0, "untested rules: {}".format(untested_rules))

    def test_parser_grammar(self):
        with open(self._test_fn, "r") as f:
            reader = csv.DictReader(f, delimiter=str("\t"))

            fail_cases = []

            for row in reader:
                if row["Func"].startswith("#"):
                    continue

                # setup input
                inputs = self._split_inputs(row["Test"], row["InType"])
                expected_results = self._split_inputs(row["Expected"],
                                                      row["InType"]) if row["Expected"] else inputs
                expected_map = dict(zip(inputs, expected_results))
                # step through each item and check
                is_valid = True if row["Valid"].lower() == "true" else False

                for key in expected_map:
                    expected_result = six.text_type(expected_map[key]).replace("u'", "'")
                    function_to_test = getattr(self.p._grammar(key), row["Func"])
                    row_str = u"{}\t{}\t{}\t{}\t{}".format(row["Func"], key, row["Valid"], "one",
                                                           expected_result)
                    try:
                        actual_result = six.text_type(function_to_test()).replace("u'", "'")
                        if not is_valid or (expected_result != actual_result):
                            print("expected: {} actual:{}".format(expected_result, actual_result))
                            fail_cases.append(row_str)
                    except Exception as e:
                        if is_valid:
                            print("expected: {} Exception: {}".format(expected_result, e))
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
            assert False, "shouldn't be here (intype = {})".format(intype)
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
