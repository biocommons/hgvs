"""Shared position handling functions for HGVS."""

import hgvs.location


def get_start_end(
    var,
) -> tuple[
    hgvs.location.SimplePosition | hgvs.location.BaseOffsetPosition,
    hgvs.location.SimplePosition | hgvs.location.BaseOffsetPosition,
]:
    """Get start and end positions from a variant or interval.

    This function handles all position types (SimplePosition, BaseOffsetPosition,
    Interval, BaseOffsetInterval) and returns the appropriate start and end positions.

    Args:
        var: A variant object with posedit.pos attribute, or an Interval object

    Returns:
        tuple: (start_position, end_position) where positions can be SimplePosition or BaseOffsetPosition
    """
    print(f"get_start_end: var: {var}")

    # Handle Interval objects directly
    if isinstance(var, hgvs.location.Interval):
        s = var.start
        e = var.end
        if not isinstance(
            s, (hgvs.location.SimplePosition, hgvs.location.BaseOffsetPosition)
        ):
            s = s.start
        if not isinstance(
            e, (hgvs.location.SimplePosition, hgvs.location.BaseOffsetPosition)
        ):
            e = e.end
        return s, e

    # Handle position objects directly
    if isinstance(
        var, (hgvs.location.SimplePosition, hgvs.location.BaseOffsetPosition)
    ):
        return var, var

    # Handle variants with posedit
    pos = var.posedit.pos

    if isinstance(pos, hgvs.location.SimplePosition):
        return pos, pos
    elif isinstance(pos, hgvs.location.BaseOffsetInterval):
        return pos.start, pos.end
    elif isinstance(pos, hgvs.location.Interval):
        s = pos.start
        e = pos.end

        if not isinstance(
            s, (hgvs.location.SimplePosition, hgvs.location.BaseOffsetPosition)
        ):
            orig_s = s
            s = s.start
            if s.is_uncertain:
                s = orig_s.end
        if not isinstance(
            e, (hgvs.location.SimplePosition, hgvs.location.BaseOffsetPosition)
        ):
            orig_e = e
            e = e.end
            if e.is_uncertain:
                e = orig_e.start
        return s, e
    else:  # BaseOffsetPosition
        return pos, pos
