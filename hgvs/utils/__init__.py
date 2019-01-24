# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function, unicode_literals

import re
from six.moves import range


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
