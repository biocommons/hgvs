import string

aa3_to_aa1_lut = {
    'Ala': 'A',    'Arg': 'R',    'Asn': 'N',    'Asp': 'D',
    'Cys': 'C',    'Gln': 'Q',    'Glu': 'E',    'Gly': 'G',
    'His': 'H',    'Ile': 'I',    'Leu': 'L',    'Lys': 'K',
    'Met': 'M',    'Phe': 'F',    'Pro': 'P',    'Ser': 'S',
    'Thr': 'T',    'Trp': 'W',    'Tyr': 'Y',    'Val': 'V',
    'Xaa': 'X',    'Ter': '*',    'Sec': 'U',
    }

aa1_to_aa3_lut = { v: k for k,v in aa3_to_aa1_lut.iteritems() }

lowercase = 'abcdefghijklmnopqrstuvwxyz'   # rather than import string



def aa3_to_aa1(s):
    "convert string of 3-letter amino acids to 1-letter amino acids"
    return None if s is None else  ''.join([aa3_to_aa1_lut[aa3] for aa3 in
                                            [ s[i:i+3] for i in range(0,len(s),3) ]])

def aa1_to_aa3(s):
    "convert string of 1-letter amino acids to 3-letter amino acids"
    return None if s is None else ''.join([aa1_to_aa3_lut[aa1] for aa1 in s])

def aa_to_aa1(s):
    "coerce string of 1- or 3-letter amino acids to 1-letter"
    return aa3_to_aa1(s) if __looks_like_aa3_p(s) else s

def aa_to_aa3(s):
    "coerce string of 1- or 3-letter amino acids to 3-letter"
    return aa1_to_aa3(s) if not __looks_like_aa3_p(s) else s


complement_transtable = string.maketrans('ACGT','TGCA')
def reverse_complement(s):
    return ''.join(reversed(s.translate(complement_transtable)))


def __looks_like_aa3_p(s):
    "string looks like a 3-letter AA string"
    return (
        s is not None
        and (len(s) % 3 == 0)
        and (len(s) == 0 or s[1] in lowercase)
        )
