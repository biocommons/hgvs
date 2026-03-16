"""Shared position handling functions for HGVS."""

import hgvs.location


def get_start_end(
    var, outer_confidence=True
) -> tuple[
    hgvs.location.SimplePosition | hgvs.location.BaseOffsetPosition,
    hgvs.location.SimplePosition | hgvs.location.BaseOffsetPosition,
]:
    """Get start and end positions from a variant or interval.

    This function handles all position types (SimplePosition, BaseOffsetPosition,
    Interval, BaseOffsetInterval) and returns the appropriate start and end positions.
    It can be expected that the returned positions have a base and an uncertain property.

    By default we return the outer confidence positions. However, if that position
    does not have a base, we return the inner confidence positions.

    TODO: add a new optional parameter that allows to define the strictness of the returned positions.
    The current behavior is more alike to an "auto" mode, since we might fall back to the inner confidence positions
    if the outer confidence positions do not have a base. A potential "strict" mode would only return the outer confidence positions,
    and raise an error if the outer confidence positions do not have a base.

    Args:
        var: A variant object with posedit.pos attribute, or an Interval object
        outer_confidence: If True, return the outer confidence positions, otherwise return the inner confidence positions

    Returns:
        tuple: (start_position, end_position) where positions can be SimplePosition or BaseOffsetPosition
    """

    # Handle Interval objects directly
    if isinstance(var, hgvs.location.Interval):
        s = var.start
        e = var.end
        if not isinstance(
            s, (hgvs.location.SimplePosition, hgvs.location.BaseOffsetPosition)
        ):
            if outer_confidence and not s.start.base is not None:
                s = s.start
            else:
                s = s.end
        if not isinstance(
            e, (hgvs.location.SimplePosition, hgvs.location.BaseOffsetPosition)
        ):
            if outer_confidence and not e.end.base is not None:
                e = e.end
            else:
                e = e.start
        return s, e

    # Handle position objects directly
    if isinstance(
        var, (hgvs.location.SimplePosition, hgvs.location.BaseOffsetPosition)
    ):
        return var, var

    # if there is no posedit, return None, all steps below this would fail
    if var.posedit is None:
        return None, None

    # Handle variants with posedit
    pos = var.posedit.pos

    if isinstance(pos, hgvs.location.SimplePosition):
        return pos, pos
    elif isinstance(pos, hgvs.location.BaseOffsetInterval):
        return pos.start, pos.end
    elif isinstance(pos, hgvs.location.Interval):
        s = pos.start
        e = pos.end

        if isinstance(s, hgvs.location.AAPosition) and isinstance(
            e, hgvs.location.AAPosition
        ):
            s = s.base
            e = e.base
            return s, e

        if not isinstance(
            s, (hgvs.location.SimplePosition, hgvs.location.BaseOffsetPosition)
        ):
            orig_s = s

            if outer_confidence and s.start.base is not None:
                s = s.start
            else:
                s = s.end
            if s.is_uncertain:
                s = orig_s.end

        if not isinstance(
            e, (hgvs.location.SimplePosition, hgvs.location.BaseOffsetPosition)
        ):
            orig_e = e
            if outer_confidence and e.end.base is not None:
                e = e.end
            else:
                e = e.start
            if e.is_uncertain:
                e = orig_e.start
        return s, e
    else:  # BaseOffsetPosition
        return pos, pos


def get_start_end_interbase(
    pos: hgvs.location.BaseOffsetPosition | hgvs.location.Interval,
    outer_confidence=True,
) -> tuple[int, int]:
    """Get start and end integer positions from a SequenceVariant.

    This function extracts integer positions from a SequenceVariant, handling uncertain positions
    and intervals. For uncertain positions, it uses the more confident boundary.

    Args:
        var: A SequenceVariant object
        outer_confidence: If True, return the outer confidence positions, otherwise return the inner confidence positions

    Returns:
        tuple: (start_int, end_int) where both are integer positions (0-based)
    """

    # Handle start position
    if pos.start.uncertain:
        # For uncertain start, use the more confident boundary
        if isinstance(pos.start, hgvs.location.Interval):
            if outer_confidence and pos.start.start.base is not None:
                seq_start = pos.start.start.base - 1
            else:
                seq_start = pos.start.end.base - 1
        else:
            seq_start = pos.start.base - 1
    else:
        # For certain start, use the start boundary
        if isinstance(pos.start, hgvs.location.Interval):
            if outer_confidence:
                seq_start = pos.start.start.base - 1
            else:
                seq_start = pos.start.end.base - 1
        else:
            seq_start = pos.start.base - 1

    # Handle end position
    if pos.end.uncertain:
        # For uncertain end, use the more confident boundary
        if isinstance(pos.end, hgvs.location.Interval):
            if outer_confidence and pos.end.end.base is not None:
                seq_end = pos.end.end.base
            else:
                seq_end = pos.end.start.base
        else:
            seq_end = pos.end.base
    else:
        # For certain end, use the end boundary
        if isinstance(pos.end, hgvs.location.Interval):
            if outer_confidence:
                seq_end = pos.end.end.base
            else:
                seq_end = pos.end.start.base
        else:
            seq_end = pos.end.base

    return seq_start, seq_end
