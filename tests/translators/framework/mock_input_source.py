#
# Mock test input source
#
from __future__ import with_statement
import csv

import hgvs.edti.interface as interface

class MockInputSource():

    def __init__(self, in_file):
        self._mock_data = self._read_input(in_file)

    def fetch_gene_info(self,ac):
        pass

    def fetch_gene_transcripts(self,ac):
        pass

    def fetch_transcript_exons(self,ac):
        return self._mock_data[ac]

    def fetch_transcript_info(self,ac):
        pass


    def get_sequence(self, accession):
        """Gets DNA sequence, given an accession; None if tag can't be converted
        """
        return self._mock_data.get(accession)

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
                                                   'cds_start': int(row['cds_start']),
                                                   'cds_stop': int(row['cds_stop']),
                                                   'protein_accession': row['protein_accession']}

        return result


def main():
    pass


if __name__ == "__main__":
    main()
