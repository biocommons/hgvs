#
# amino acid translation utilities
#

CODES_1_3 = {'A': 'Ala',
             'R': 'Arg',
             'N': 'Asn',
             'D': 'Asp',
             'C': 'Cys',
             'E': 'Glu',
             'Q': 'Gln',
             'G': 'Gly',
             'H': 'His',
             'I': 'Ile',
             'L': 'Leu',
             'K': 'Lys',
             'M': 'Met',
             'F': 'Phe',
             'P': 'Pro',
             'S': 'Ser',
             'T': 'Thr',
             'W': 'Trp',
             'Y': 'Tyr',
             'V': 'Val',
             'X': 'Xaa',
             '*': 'Ter'
             }


def convert_AA_1_to_3(aa_seq):
    """convert single char AA notation to 3 char; also converts * and X (to Ter and Xaa)
    :param aa_seq: string of amino acids represented by 1 char
    :type iterable
    :return amino acids string represented by 3 char codes"""
    result = ''.join(CODES_1_3[x] for x in aa_seq)
    return result


def main():
    pass


if __name__ == "__main__":
    main()

