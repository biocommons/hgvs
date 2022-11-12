import pytest

from hgvs.exceptions import HGVSInvalidIntervalError
from hgvs.utils.cigarmapper import CIGARMapper

cigar = "3=2N=X=3N=I=D="
#  0   1   2           3   4   5               6       7   8   9  tgt
#  =   =   =   N   N   =   X   =   N   N   N   =   I   =   D   =
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13      14  ref


cm = CIGARMapper(cigar)


def test_cigarmapper():
    """simple cigar test"""

    assert cm.ref_len == 15
    assert cm.tgt_len == 10

    assert len(cm.ref_pos) == len(cm.tgt_pos)
    assert cm.ref_pos == [0, 3, 5, 6, 7, 8, 11, 12, 13, 14, 14, 15]
    assert cm.tgt_pos == [0, 3, 3, 4, 5, 6, 6, 7, 7, 8, 9, 10]

    # ref to tgt
    assert cm.map_ref_to_tgt(0, "start") == (0, 0, "=")
    assert cm.map_ref_to_tgt(0, "end") == (0, 0, "=")
    assert cm.map_ref_to_tgt(1, "start") == (1, 0, "=")
    assert cm.map_ref_to_tgt(1, "end") == (1, 0, "=")
    assert cm.map_ref_to_tgt(2, "start") == (2, 0, "=")
    assert cm.map_ref_to_tgt(2, "end") == (2, 0, "=")
    assert cm.map_ref_to_tgt(3, "start") == (2, 1, "N")
    assert cm.map_ref_to_tgt(3, "end") == (2, 1, "N")
    assert cm.map_ref_to_tgt(4, "start") == (3, -1, "N")
    assert cm.map_ref_to_tgt(4, "end") == (3, -1, "N")
    assert cm.map_ref_to_tgt(5, "start") == (3, 0, "=")
    assert cm.map_ref_to_tgt(5, "end") == (3, 0, "=")
    assert cm.map_ref_to_tgt(6, "start") == (4, 0, "X")
    assert cm.map_ref_to_tgt(6, "end") == (4, 0, "X")
    assert cm.map_ref_to_tgt(7, "start") == (5, 0, "=")
    assert cm.map_ref_to_tgt(7, "end") == (5, 0, "=")
    assert cm.map_ref_to_tgt(8, "start") == (5, 1, "N")
    assert cm.map_ref_to_tgt(8, "end") == (5, 1, "N")
    assert cm.map_ref_to_tgt(9, "start") == (5, 2, "N")
    assert cm.map_ref_to_tgt(9, "end") == (5, 2, "N")
    assert cm.map_ref_to_tgt(10, "start") == (6, -1, "N")
    assert cm.map_ref_to_tgt(10, "end") == (6, -1, "N")
    assert cm.map_ref_to_tgt(11, "start") == (6, 0, "=")
    assert cm.map_ref_to_tgt(11, "end") == (6, 0, "=")
    assert cm.map_ref_to_tgt(12, "start") == (6, 0, "I")
    assert cm.map_ref_to_tgt(12, "end") == (7, 0, "I")
    assert cm.map_ref_to_tgt(13, "start") == (7, 0, "=")
    assert cm.map_ref_to_tgt(13, "end") == (7, 0, "=")
    assert cm.map_ref_to_tgt(14, "start") == (9, 0, "=")
    assert cm.map_ref_to_tgt(14, "end") == (9, 0, "=")

    # tgt to ref
    assert cm.map_tgt_to_ref(0, "start") == (0, 0, "=")
    assert cm.map_tgt_to_ref(0, "end") == (0, 0, "=")
    assert cm.map_tgt_to_ref(1, "start") == (1, 0, "=")
    assert cm.map_tgt_to_ref(1, "end") == (1, 0, "=")
    assert cm.map_tgt_to_ref(2, "start") == (2, 0, "=")
    assert cm.map_tgt_to_ref(2, "end") == (2, 0, "=")
    assert cm.map_tgt_to_ref(3, "start") == (5, 0, "=")
    assert cm.map_tgt_to_ref(3, "end") == (5, 0, "=")
    assert cm.map_tgt_to_ref(4, "start") == (6, 0, "X")
    assert cm.map_tgt_to_ref(4, "end") == (6, 0, "X")
    assert cm.map_tgt_to_ref(5, "start") == (7, 0, "=")
    assert cm.map_tgt_to_ref(5, "end") == (7, 0, "=")
    assert cm.map_tgt_to_ref(6, "start") == (11, 0, "=")
    assert cm.map_tgt_to_ref(6, "end") == (11, 0, "=")
    assert cm.map_tgt_to_ref(7, "start") == (13, 0, "=")
    assert cm.map_tgt_to_ref(7, "end") == (13, 0, "=")
    assert cm.map_tgt_to_ref(8, "start") == (13, 0, "D")
    assert cm.map_tgt_to_ref(8, "end") == (14, 0, "D")
    assert cm.map_tgt_to_ref(9, "start") == (14, 0, "=")
    assert cm.map_tgt_to_ref(9, "end") == (14, 0, "=")


def test_cigarmapper_strict_bounds():
    # exception raised for out of bounds on left?
    with pytest.raises(HGVSInvalidIntervalError):
        cm.map_ref_to_tgt(-1, "start", strict_bounds=True)

    # ... and right?
    with pytest.raises(HGVSInvalidIntervalError):
        cm.map_ref_to_tgt(cm.ref_len + 1, "start", strict_bounds=True)

    # test whether 1 base outside bounds results in correct position
    assert (0, 0, "=") == cm.map_ref_to_tgt(0, "start", strict_bounds=True)
    assert (-1, 0, "=") == cm.map_ref_to_tgt(-1, "start", strict_bounds=False)
    assert (cm.tgt_len, 0, "=") == cm.map_ref_to_tgt(cm.ref_len, "start", strict_bounds=True)
    assert (cm.tgt_len - 1, 0, "=") == cm.map_ref_to_tgt(cm.ref_len - 1, "start", strict_bounds=False)
