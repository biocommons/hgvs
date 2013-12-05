"""
stopgap -- components that need to be fleshed out later

importing stopgap should tell you that you're using incomplete and/or unstable code
"""

import hashlib

def pseq_to_ac(s):
    "placeholder function that will return an accession for a given protein sequence"
    return 'MD5_' + hashlib.md5(s).hexdigest()[:8]
