import base64
import hashlib
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


def seq_seguid(seq):
    """returns BioPython-compatible seguid""" 
    seq = seq.upper().rstrip('*')
    return base64.standard_b64encode(hashlib.sha1(seq).digest()).rstrip('=')

def seq_md5(seq):
    "returns sequence md5 as hex digest"
    seq = seq.upper().rstrip('*')
    return hashlib.md5(seq).hexdigest()

def seq_sha1(seq):
    """returns sequence sha1 as url-safe base64 encoded digest. This is
    identical to seguid except for hashes that contain / or +."""
    seq = seq.upper().rstrip('*')
    return base64.urlsafe_b64encode(hashlib.sha1(seq).digest()).rstrip('=')


def __looks_like_aa3_p(s):
    "string looks like a 3-letter AA string"
    return (
        s is not None
        and (len(s) % 3 == 0)
        and (len(s) == 0 or s[1] in string.lowercase)
        )

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
