# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

# Mock test input source

import unicodecsv as csv


class MockInputSource():
    def __init__(self, in_file):
        self._mock_data = self._read_input(in_file)

    def fetch_gene_info(self, ac):
        pass

    def fetch_gene_transcripts(self, ac):
        pass

    def fetch_transcript_exons(self, ac):
        result = None
        data = self._mock_data.get(ac)
        if data:
            result = {
                'ord': 1,
                't_start_i': 0,
                't_end_i': data['cds_end_i'] - data['cds_start_i'],
                't_seq_a': data['transcript_sequence']
            }

        return [result]

    def fetch_transcript_info(self, ac):
        result = None
        data = self._mock_data.get(ac)
        if data:    # interbase coordinates
            result = {'cds_start_i': data['cds_start_i'], 'cds_end_i': data['cds_end_i']}
        return result

    def get_tx_identity_info(self, ac):
        return self.fetch_transcript_info(ac)

    def get_tx_info(self, ac):
        return self.fetch_transcript_info(ac)

    def get_tx_exons(self, ac):
        return self.fetch_transcript_exons(ac)

    def get_tx_seq(self, ac):
        result = None
        data = self._mock_data.get(ac)
        if data:
            result = data['transcript_sequence']
        return result

    def fetch_seq(self, ac, start_i=None, end_i=None):
        return self.get_tx_seq(ac)[start_i:end_i]

    #
    # internal methods
    #

    def _read_input(self, in_file):
        """Dummy file of inputs

        :param in_file: path to input file of 2 cols (tab-delim); accession_number, sequence
        :type string
        :return dictionary of accession_number to sequence tags
        """
        result = {}
        with open(in_file, 'r') as f:
            reader = csv.DictReader(f, delimiter=str('\t'))
            for row in reader:
                result[row['accession']] = {
                    'transcript_sequence': row['transcript_sequence'],
                    'cds_start_i': int(row['cds_start_i']),
                    'cds_end_i': int(row['cds_end_i'])
                }

        return result


def main():
    pass


if __name__ == "__main__":
    main()

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
