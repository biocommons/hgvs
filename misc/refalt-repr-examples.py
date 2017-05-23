#!/usr/bin/env python

import hgvs
import hgvs.parser
from tabulate import tabulate
from six.moves import map

hp = hgvs.parser.Parser()

variants = [
    "NM_01.2:c.1_2del",
    "NM_01.2:c.1_2del2",
    "NM_01.2:c.1_2delAA",
    "NM_01.2:c.1_2insTTT",
    "NM_01.2:c.1_2delinsTTT",
    "NM_01.2:c.1_2del2insTTT",
    "NM_01.2:c.1_2delAAinsTTT",
    "NM_01.2:c.1_2delAAinsAA",
    "NM_01.2:c.1_1delAinsA",
    "NM_01.2:c.1A>A",
    "NM_01.2:c.1_1delAinsT",
    "NM_01.2:c.1A>T",
    "NM_01.2:c.1=",
    "NM_01.2:c.1A=",
    "NM_01.2:c.1AA=",
]


headers = "hgvs ref alt etype delta".split()
def gen1(h):
    v = hp.parse_hgvs_variant(h)
    return [h, v.posedit.edit.ref, v.posedit.edit.alt, v.posedit.edit.type, v.posedit.length_change()]


rows = [list(map(str, gen1(h))) for h in variants]

print(tabulate(rows, headers=headers)) #, tablefmt="pipe"))
