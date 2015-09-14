#!/usr/bin/env python
"""demonstration script for bug #274
Thanks to Peter Causey-Freeman

Actually, two bugs:

1) Using n_to_c (or g_to_c) with a non-coding transcript resulted in a
type error. An explicit warning is now generated.

2) In EVM, the normalizer was always initialized with the default
alt_aln_method.

"""

import hgvs
import hgvs.dataproviders.uta 
import hgvs.variantmapper
import hgvs.parser 
from hgvs.exceptions import HGVSUsageError

hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()
alt_aln_method = "genebuild"
primary_assembly = 'GRCh37'
variant = 'ENST00000225964:c.589G>T'
evm = hgvs.variantmapper.EasyVariantMapper(hdp, primary_assembly, alt_aln_method)
hgvs_variant = hp.parse_hgvs_variant(variant)
hgvs_genomic = evm.c_to_g(hgvs_variant)
rts = evm.relevant_transcripts(hgvs_genomic)


tx_ac = 'ENST00000495677'
try:
    v1 = evm.g_to_c(hgvs_genomic, tx_ac)
    print(v1)
except HGVSUsageError as e:
    print(e)

v1b = evm.g_to_n(hgvs_genomic, tx_ac)
print(v1b)

tx_ac = 'ENST00000225964'
v2 = evm.g_to_c(hgvs_genomic, tx_ac)
print(v2)
