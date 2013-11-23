aa3_to_aa1_lut = {
    'Ala': 'A',    'Arg': 'R',    'Asn': 'N',    'Asp': 'D',
    'Cys': 'C',    'Gln': 'Q',    'Glu': 'E',    'Gly': 'G',
    'His': 'H',    'Ile': 'I',    'Leu': 'L',    'Lys': 'K',
    'Met': 'M',    'Phe': 'F',    'Pro': 'P',    'Ser': 'S',
    'Thr': 'T',    'Trp': 'W',    'Tyr': 'Y',    'Val': 'V'
    }

aa1_to_aa3_lut = { v: k for k,v in aa3_to_aa1_lut.iteritems() }

lowercase = 'abcdefghijklmnopqrstuvwxyz'   # rather than import string



def aa3_to_aa1(s):
    "convert string of 3-letter amino acids to 1-letter amino acids"
    return ''.join([aa3_to_aa1_lut[aa3] for aa3 in
                    [ s[i:i+3] for i in range(0,len(s),3) ]])

def aa1_to_aa3(s):
    "convert string of 1-letter amino acids to 3-letter amino acids"
    return ''.join([aa1_to_aa3_lut[aa1] for aa1 in s])

def aa_to_aa1(s):
    "coerce string of 1- or 3-letter amino acids to 1-letter"
    if s[1] in lowercase:
        return aa3_to_aa1(s)
    return s

def aa_to_aa3(s):
    "coerce string of 1- or 3-letter amino acids to 3-letter"
    if s[1] in lowercase:
        return s
    return aa1_to_aa3(s)
