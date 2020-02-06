# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function, unicode_literals

import itertools
import re


def build_tx_cigar(exons, strand):
    """builds a single CIGAR string representing an alignment of the
    transcript sequence to a reference sequence, including introns.
    The input exons are expected to be in transcript order, and the
    resulting CIGAR is also in transcript order. 

    >>> build_tx_cigar([], 1) is None
    True
    """
    cigarelem_re = re.compile(r"\d+[=DIMNX]")

    def _reverse_cigar(c):
        return ''.join(reversed(cigarelem_re.findall(c)))

    if len(exons) == 0:
        return None

    # flip orientation of all CIGARs if on - strand
    if strand == -1:
        cigars = [_reverse_cigar(e["cigar"]) for e in exons]
    else:
        cigars = [e["cigar"] for e in exons]

    tx_cigar = [cigars[0]]    # exon 1
    for i in range(1, len(cigars)):    # and intron + exon pairs thereafter
        intron = str(exons[i]["alt_start_i"] - exons[i - 1]["alt_end_i"]) + "N"
        tx_cigar += [intron, cigars[i]]

    tx_cigar_str = "".join(tx_cigar)

    return tx_cigar_str



def parse_cigar(cigar):
    """For a given CIGAR string, return the start positions of each
    aligned segment in ref and tgt, and a list of CIGAR operators.

    """
    cigar_re = re.compile(r"(?P<len>\d+)(?P<op>[=DIMNX])")
    ces = [m.groupdict() for m in cigar_re.finditer(cigar)]
    ref_pos = [None] * len(ces)
    tgt_pos = [None] * len(ces)
    cigar_op = [None] * len(ces)
    ref_cur = tgt_cur = 0
    for i, ce in enumerate(ces):
        ref_pos[i] = ref_cur
        tgt_pos[i] = tgt_cur
        cigar_op[i] = ce["op"]
        step = int(ce["len"])
        if ce["op"] in "=MINX":
            ref_cur += step
        if ce["op"] in "=MDX":
            tgt_cur += step
    ref_pos.append(ref_cur)
    tgt_pos.append(tgt_cur)
    return ref_pos, tgt_pos, cigar_op




# <LICENSE>
# Copyright 2018 HGVS Contributors (https://github.com/biocommons/hgvs)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
