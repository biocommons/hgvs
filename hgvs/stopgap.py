"""
stopgap -- components that need to be fleshed out later

importing stopgap should tell you that you're using incomplete and/or unstable code
"""

import gzip
import json
from pkg_resources import resource_filename

from Bio.SeqUtils.CheckSum import seguid

seguid_acs = json.loads(gzip.open(resource_filename(__name__,'/data/seguid-acs.json.gz')).read())

def pseq_to_ac(seq):
    "placeholder function that returns an accession for a given protein sequence"
    seq = seq.rstrip('*')
    hash = seguid(seq)
    try:
        return seguid_acs[hash][0]
    except KeyError:
        return 'SEGUID_' + hash
