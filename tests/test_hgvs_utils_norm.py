import random

import pytest

from hgvs.utils.norm import normalize_alleles_left, normalize_alleles_right

#        01234567
REF = "TCACACAG"

SHUFFLES = [
    # (normalize function, start, stop, alleles, bound, ref_step, expected)
    (normalize_alleles_left, 7, 7, ("", "CA"), 0, 2, (1, 1, ("", "CA"))),
    (normalize_alleles_left, 6, 6, ("", "AC"), 0, 2, (1, 1, ("", "CA"))),  # rotated
    (normalize_alleles_left, 6, 6, ("", "CAC"), 0, 2, (3, 3, ("", "CAC"))),
    (normalize_alleles_left, 6, 6, ("", "ACA"), 0, 2, (6, 6, ("", "ACA"))),  # no shift
    (normalize_alleles_left, 5, 7, ("CA", ""), 0, 3, (1, 3, ("CA", ""))),
    (normalize_alleles_left, 7, 7, ("", "CA"), 3, 2, (3, 3, ("", "CA"))),  # stops at bound
    (normalize_alleles_right, 1, 1, ("", "CA"), 8, 2, (7, 7, ("", "CA"))),
    (normalize_alleles_right, 2, 2, ("", "ACA"), 8, 2, (5, 5, ("", "ACA"))),
    (normalize_alleles_right, 1, 3, ("CA", ""), 8, 3, (5, 7, ("CA", ""))),
]


@pytest.mark.parametrize("shuffle", SHUFFLES)
def test_shuffle_indel(shuffle):
    func, *args, expected = shuffle
    assert tuple(func(REF, *args)) == expected


def test_shuffle_large_dup():
    """Shuffling an insertion past a 1Mb copy of itself was quadratic (minutes per Mb)"""
    rng = random.Random(0)  # noqa: S311
    dup = "".join(rng.choice("ACGT") for _ in range(1_000_000))
    # Padding that can't shuffle any further
    left_pad = "A" if dup[-1] != "A" else "C"
    right_pad = "A" if dup[0] != "A" else "C"
    ref = left_pad + dup + right_pad
    end = len(dup) + 1
    left = normalize_alleles_left(ref, end, end, ("", dup), 0, 20)
    assert tuple(left) == (1, 1, ("", dup))
    right = normalize_alleles_right(ref, 1, 1, ("", dup), len(ref), 20)
    assert tuple(right) == (end, end, ("", dup))
