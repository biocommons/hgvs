#
# Mock test input source
#
from __future__ import with_statement
import csv

import hgvs.edti.interface as interface

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
            result = {'ord': 1,
                      't_start_i': 0,
                      't_end_i': data['cds_stop_i'] - data['cds_start_i'],
                      't_seq_a': data['transcript_sequence']
            }

        return [result]

    def fetch_transcript_info(self,ac):
        result = None
        data = self._mock_data.get(ac)
        if data:     # interbase coordinates
            result = {'cds_start_i': data['cds_start_i'],
                      'cds_stop_i': data['cds_stop_i']}
        return result

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
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                result[row['accession']] = {'transcript_sequence': row['transcript_sequence'],
                                                   'cds_start_i': int(row['cds_start_i']),
                                                   'cds_stop_i': int(row['cds_stop_i'])}

        return result


def main():
    pass


if __name__ == "__main__":
    main()
