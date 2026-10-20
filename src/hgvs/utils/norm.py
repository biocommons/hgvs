"""Utility to normalize variants.
Code of the normalization utilities were imported from vgraph
https://github.com/bioinformed/vgraph

"""

from collections import namedtuple

from bioutils.normalize import roll_left, roll_right


def trim_common_suffixes(strs, min_len=0):
    """
    trim common suffixes

    >>> trim_common_suffixes("A", 1)
    (0, 'A')

    """

    if len(strs) < 2:
        return 0, strs

    rev_strs = [s[::-1] for s in strs]

    trimmed, rev_strs = trim_common_prefixes(rev_strs, min_len)

    if trimmed:
        strs = [s[::-1] for s in rev_strs]

    return trimmed, strs


def trim_common_prefixes(strs, min_len=0):
    """trim common prefixes"""

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


def _rotate(seq, k):
    """Rotate seq left by k bases (k may be negative or exceed len(seq))"""
    if not seq:
        return seq
    k %= len(seq)
    return seq[k:] + seq[:k]


def normalize_alleles_left(ref, start, stop, alleles, bound, ref_step, shuffle=True):  # noqa: ARG001
    """
    Normalize loci by removing extraneous reference padding

    ref_step is unused and kept for compatibility

    >>> normalize_alleles_left("A", 1, 2, "A", 1, 2)
    shuffled_alleles(start=1, stop=2, alleles='A')

    """

    normalized_alleles = namedtuple("shuffled_alleles", "start stop alleles")

    if len(alleles) < 2:
        return normalized_alleles(start, stop, alleles)

    # STEP 1: Trim common suffix
    trimmed, alleles = trim_common_suffixes(alleles)
    stop -= trimmed

    # STEP 2: Trim common prefix
    trimmed, alleles = trim_common_prefixes(alleles)
    start += trimmed

    # assert bound <= start,'start={:d}, left bound={:d}'.format(start, bound)

    # STEP 3: Left shuffle. Shifting an insertion or deletion one base rotates its sequence,
    #         so bioutils counts the shift with a cyclic comparison and the alleles are
    #         rotated once (rebuilding them every ref_step bases was quadratic for a large dup)
    if shuffle:
        shift = roll_left(ref.upper(), alleles, start - 1, bound)
        if shift:
            alleles = [_rotate(a, -shift) for a in alleles]
            start -= shift
            stop -= shift

    return normalized_alleles(start, stop, tuple(alleles))


def normalize_alleles_right(ref, start, stop, alleles, bound, ref_step, shuffle=True):  # noqa: ARG001
    """
    Normalize loci by removing extraneous reference padding

    ref_step is unused and kept for compatibility

    >>> normalize_alleles_right("A", 1, 2, "A", 1, 2)
    shuffled_alleles(start=1, stop=2, alleles='A')

    """

    normalized_alleles = namedtuple("shuffled_alleles", "start stop alleles")

    if len(alleles) < 2:
        return normalized_alleles(start, stop, alleles)

    # STEP 1: Trim common prefix
    trimmed, alleles = trim_common_prefixes(alleles)
    start += trimmed

    # STEP 2: Trim common suffix
    trimmed, alleles = trim_common_suffixes(alleles)
    stop -= trimmed

    # assert bound >= stop,'stop={:d}, right bound={:d}'.format(stop, bound)

    # STEP 3: Right shuffle, see normalize_alleles_left
    if shuffle:
        shift = roll_right(ref.upper(), alleles, stop, bound - 1)
        if shift:
            alleles = [_rotate(a, shift) for a in alleles]
            start += shift
            stop += shift

    return normalized_alleles(start, stop, tuple(alleles))


def normalize_alleles(ref, start, stop, alleles, bound, ref_step, left, shuffle=True):
    if left:
        return normalize_alleles_left(ref, start, stop, alleles, bound, ref_step, shuffle)
    return normalize_alleles_right(ref, start, stop, alleles, bound, ref_step, shuffle)


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
