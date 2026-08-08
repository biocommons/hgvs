"""Regression tests for a normalizer bug uncovered while fixing the
AttributeError crash in hgvs.utils.position.get_start_end() (a follow-up
to #801).

#801 taught the c./n. grammar to parse uncertain intervals like
c.(?_6)_245. Every side of one of these intervals is built via
def_c_interval/def_n_interval, which *always* wraps its result in a
BaseOffsetInterval -- even a "plain" position like `245` becomes a
degenerate BaseOffsetInterval(start=245, end=245, uncertain=False)
rather than a bare BaseOffsetPosition.

hgvs.normalizer.Normalizer.normalize() (normalizer.py, around L210-214)
shuffles a variant's boundary by mutating a single position object in
place:

    s_norm, e_norm = get_start_end(var_norm)
    ...
    s_norm.base = ref_start
    e_norm.base = ref_end

get_start_end() correctly unwraps the nested interval and returns one of
its two aliased positions, but only that one gets the shuffled value.
The other half of the degenerate pair is left stale, so the two fall out
of sync: the printed HGVS interval gains a spurious extra number (e.g.
"(?_7)_245_246del" instead of "(?_7)_246del"), or, when both start and
end are affected, the resulting interval fails its own start<=end
sanity check and normalize() raises instead.

These tests are marked xfail(strict=True) until the normalizer stops
mutating a single aliased position in place: they should fail today,
and start failing loudly (XPASS) the moment someone fixes it, as a
prompt to remove the marker. See PR discussion for #801.
"""

import pytest

import hgvs.normalizer


@pytest.fixture(scope="module")
def norm(hdp):
    return hgvs.normalizer.Normalizer(hdp, shuffle_direction=3, cross_boundaries=True)


# Each of these wraps a position that is *known* to shift under 3'
# normalization when written as a plain, non-uncertain deletion, e.g.
# NM_000333.3:c.245del -> c.246del and NM_001166478.1:c.31del -> c.35del
# (both already covered by test_hgvs_normalizer.py, unaffected by this
# bug). Writing the same position as one side of an uncertain interval
# is what triggers the desync.
degenerate_sync_cases = [
    "NM_000333.3:c.(?_6)_245del",
    "NM_001166478.1:c.(?_31)_31del",
]


@pytest.mark.xfail(
    reason="normalizer.py mutates one aliased half of a degenerate nested "
    "position in place, leaving the other stale (follow-up to #801)",
    strict=True,
)
@pytest.mark.parametrize("hgvs_str", degenerate_sync_cases)
def test_normalize_keeps_degenerate_nested_position_in_sync(hgvs_str, parser, norm):
    """A plain (non-uncertain) side of an uncertain c./n. interval is
    internally a degenerate BaseOffsetInterval(pos, pos). Normalizing
    must keep both halves of that pair equal -- otherwise the printed
    HGVS string gains a spurious third number.
    """
    var = parser.parse(hgvs_str)
    var_norm = norm.normalize(var)

    for side in (var_norm.posedit.pos.start, var_norm.posedit.pos.end):
        if hasattr(side, "start") and hasattr(side, "end") and not side.uncertain:
            assert side.start.base == side.end.base, (
                f"{hgvs_str}: degenerate nested position fell out of sync "
                f"after normalization -> {var_norm}"
            )


# These hit the same root cause but end up with a start/end ordering
# that fails Interval.validate()'s own start<=end sanity check, so
# normalize() raises HGVSInvalidVariantError outright instead of
# printing a malformed-but-parseable string.
crash_cases = [
    "NM_001166478.1:c.31_(31_31)del",
    "NM_000088.3:c.(?_589)_600inv",
]


@pytest.mark.xfail(
    reason="normalizer.py mutates one aliased half of a degenerate nested "
    "position in place, leaving the other stale (follow-up to #801)",
    strict=True,
)
@pytest.mark.parametrize("hgvs_str", crash_cases)
def test_normalize_does_not_raise_on_uncertain_interval(hgvs_str, parser, norm):
    """Normalizing a valid, parseable uncertain c./n. interval must not
    raise.
    """
    var = parser.parse(hgvs_str)
    norm.normalize(var)  # should not raise
