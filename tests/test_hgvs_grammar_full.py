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
# - list: input is a list delimited by a pipe character ('|')
# Expected: expected result (if stringifying input does not return the same answer, e,g. '+1' -> '1')
# - if expected is left blank, then it is assumed that stringifying the parsed input returns the same answer.
#
from __future__ import with_statement

import os
import csv
import pprint
import unittest

import hgvs.parser

class TestGrammarFull(unittest.TestCase):

    def setUp(self):
        self.p = hgvs.parser.Parser()
        self.grammar = self.p._grammar

    def test_parser_grammar(self):
        fn = os.path.join(os.path.dirname(__file__), 'data', 'grammar_test.tsv')
        with open(fn, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')

            fail_cases = []

            for row in reader:
                if row['Func'].startswith('#'):
                    continue

                # setup input
                inputs = self._split_inputs(row['Test'], row['InType'])
                expected_results = self._split_inputs(row['Expected'], row['InType']) \
                    if row['Expected'] else inputs
                expected_map = dict(zip(inputs, expected_results))

                # step through each item and check
                is_valid = True if row['Valid'].lower() == 'true' else False

                for key in expected_map:
                    expected_result = str(expected_map[key])
                    function_to_test = getattr(self.p._grammar(key), row['Func'])
                    row_str = "{}\t{}\t{}\t{}\t{}".format(row['Func'], key, row['Valid'], 'one', expected_result)
                    try:
                        actual_result = str(function_to_test())
                        if not is_valid or (expected_result != actual_result):
                            fail_cases.append(row_str)
                    except Exception as e:
                        if is_valid:
                            fail_cases.append(row_str)

        # everything should have passed - report whatever failed
        self.assertTrue(len(fail_cases) == 0, pprint.pprint(fail_cases))


    def _split_inputs(self, in_string, intype):
        DELIM = '|'
        if intype == 'list':
            inputs = in_string.split(DELIM)
        elif intype == 'string':
            inputs = list(in_string)
        else:   # intype == 'one'
            inputs = [in_string]
        inputs = [x if x != 'None' else None for x in inputs]
        return inputs

if __name__ == '__main__':
    unittest.main()

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
