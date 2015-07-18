# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

#
# Utility to normalize variants.
# Code of the normalization utilities were imported from vgraph
# https://github.com/bioinformed/vgraph
#

from collections import namedtuple


def trim_common_suffixes(strs, min_len=0):
    '''trim common suffixes'''

    if len(strs) < 2:
        return 0, strs

    rev_strs = [s[::-1] for s in strs]

    trimmed, rev_strs = trim_common_prefixes(rev_strs, min_len)

    if trimmed:
        strs = [s[::-1] for s in rev_strs]

    return trimmed, strs


def trim_common_prefixes(strs, min_len=0):
    '''trim common prefixes'''

    trimmed = 0

    if len(strs) > 1:
        s1 = min(strs)
        s2 = max(strs)

        for i in range(len(s1) - min_len):
            if s1[i] != s2[i]:
                break
            trimmed = i + 1

    if trimmed > 0:
        strs = [s[trimmed:] for s in strs]

    return trimmed, strs


def normalize_alleles_left(ref, start, stop, alleles, bound, ref_step, shuffle=True):
    '''Normalize loci by removing extraneous reference padding'''

    normalized_alleles = namedtuple('shuffled_alleles', 'start stop alleles')

    if len(alleles) < 2 or start <= 0 or stop <= 0:
        return normalized_alleles(start, stop, alleles)

    # STEP 1: Trim common suffix
    trimmed, alleles = trim_common_suffixes(alleles)
    stop -= trimmed

    # STEP 2: Trim common prefix
    trimmed, alleles = trim_common_prefixes(alleles)
    start += trimmed

    #assert bound <= start,'start={:d}, left bound={:d}'.format(start, bound)

    # STEP 3: While a null allele exists, left shuffle by prepending alleles
    #         with reference and trimming common suffixes
    while shuffle and '' in alleles and start > bound:
        step = min(ref_step, start - bound)

        r = ref[start - step:start].upper()
        new_alleles = [r + a for a in alleles]

        trimmed, new_alleles = trim_common_suffixes(new_alleles)

        if not trimmed:
            break

        start -= trimmed
        stop -= trimmed

        if trimmed == step:
            alleles = new_alleles
        else:
            left = step - trimmed
            alleles = [a[left:] for a in new_alleles]
            break

    return normalized_alleles(start, stop, tuple(alleles))


def normalize_alleles_right(ref, start, stop, alleles, bound, ref_step, shuffle=True):
    '''Normalize loci by removing extraneous reference padding'''

    normalized_alleles = namedtuple('shuffled_alleles', 'start stop alleles')

    chrom_stop = len(ref)

    if len(alleles) < 2 or stop >= chrom_stop:
        return normalized_alleles(start, stop, alleles)

    # STEP 1: Trim common prefix
    trimmed, alleles = trim_common_prefixes(alleles)
    start += trimmed

    # STEP 2: Trim common suffix
    trimmed, alleles = trim_common_suffixes(alleles)
    stop -= trimmed

    #assert bound >= stop,'stop={:d}, right bound={:d}'.format(stop, bound)

    # STEP 3: While a null allele exists, right shuffle by appending alleles
    #         with reference and trimming common prefixes
    while shuffle and '' in alleles and stop < bound:
        step = min(ref_step, bound - stop)

        r = ref[stop:stop + step].upper()
        new_alleles = [a + r for a in alleles]

        trimmed, new_alleles = trim_common_prefixes(new_alleles)

        if not trimmed:
            break

        start += trimmed
        stop += trimmed

        if trimmed == step:
            alleles = new_alleles
        else:
            left = step - trimmed
            alleles = [a[:-left] for a in new_alleles]
            break

    return normalized_alleles(start, stop, tuple(alleles))


def normalize_alleles(ref, start, stop, alleles, bound, ref_step, left, shuffle=True):
    if left:
        return normalize_alleles_left(ref, start, stop, alleles, bound, ref_step, shuffle)
    else:
        return normalize_alleles_right(ref, start, stop, alleles, bound, ref_step, shuffle)

## <LICENSE>
## Copyright 2015 Kevin Jacobs (https://github.com/bioinformed/vgraph)
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
