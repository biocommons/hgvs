import string

NC_to_chr_dict = {
    'NC_000001.10':  '1', 'NC_000002.11':  '2', 'NC_000003.11': '3',
    'NC_000004.11':  '4', 'NC_000005.9' :  '5', 'NC_000006.11': '6',
    'NC_000007.13':  '7', 'NC_000008.10':  '8', 'NC_000009.11': '9',
    'NC_000010.10': '10', 'NC_000011.9' : '11', 'NC_000012.11': '12',
    'NC_000013.10': '13', 'NC_000014.8' : '14', 'NC_000015.9' : '15',
    'NC_000016.9' : '16', 'NC_000017.10': '17', 'NC_000018.9' : '18',
    'NC_000019.9' : '19', 'NC_000020.10': '20', 'NC_000021.8' : '21',
    'NC_000022.10': '22', 'NC_000023.10':  'X', 'NC_000024.9' : 'Y',
    }

chr_to_NC_dict = dict([ (v,k) for k,v in NC_to_chr_dict.iteritems() ])

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
    return ''.join(reversed(s.translate(complement_transtable))) if s is not None else None


def __looks_like_aa3_p(s):
    "string looks like a 3-letter AA string"
    return (
        s is not None
        and (len(s) % 3 == 0)
        and (len(s) == 0 or s[1] in lowercase)
        )

def chr_to_nc(s):
    return chr_to_NC_dict[s]
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
