#!/usr/bin/env python

import hgvs
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.sequencevariant
import hgvs.assemblymapper

hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()
evm = hgvs.assemblymapper.AssemblyMapper(hdp,
                                           replace_reference=True, assembly_name='GRCh37',
                                           alt_aln_method='splign')

#v = hp.parse_hgvs_variant("NM_000059.3:c.7790A>G")
#v = hp.parse_hgvs_variant("NM_000059.3:c.7790_7792delAAG")
v = hp.parse_hgvs_variant("NM_000059.3:c.7790delAAG")
evm.c_to_p(v)

