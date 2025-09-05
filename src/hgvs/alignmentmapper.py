# -*- coding: utf-8 -*-
"""Mapping positions between pairs of sequence alignments

The AlignmentMapper class is at the heart of mapping between aligned sequences.

"""

# Implementation note re: "no-zero correction": HGVS doesn't have a
# 0. Counting is -3, -2, -1, 1, 2, 3 :-/ Coordinate calculations must
# take this discontinuity in c. positions into account.  The
# implementation of imaginary transcript positions creates a second
# discontinuity. (By analogy with c., n.0 is declared to not exist.)
# The strategy used in this code is to use internal c0 and n0
# coordinates, which include 0, for coordinate calculations and to
# translate these to c. and n. positions as needed.
#
#                imag.                                                 imag.
#              upstream     5' UTR           CDS          3' UTR      downstr
#                                     |>
#            - - - - - - ———————————— ||||||||||||||||| ——————————— - - - - - -
#                           a     b     C     D     E     f     g     h     i
#    c.        -4    -3    -2    -1  !  1     2     3  ! *1    *2    *3    *4
#    c0        -4    -3    -2    -1     0     1     2     3     4     5     6
#    n0        -2    -1     0     1     2     3     4     5     6     7     8
#    n.        -2    -1  !  1     2     3     4     5     6     7     8     9
#    g.   ... 123   124   125   126   127   128   129   130   131   132   133 ...
#

from __future__ import absolute_import, division, print_function, unicode_literals

from typing import Optional

from bioutils.coordinates import strand_int_to_pm

import hgvs.location
from hgvs import global_config
from hgvs.enums import Datum
from hgvs.exceptions import (
    HGVSDataNotAvailableError,
    HGVSError,
    HGVSInvalidIntervalError,
    HGVSUsageError,
)
from hgvs.location import Interval, BaseOffsetInterval
from hgvs.utils import build_tx_cigar
from hgvs.utils.cigarmapper import CIGARMapper
import math


def _zbc_to_hgvs(i: int):
    """Convert zero-based coordinate to hgvs (1 based, missing zero)"""
    if i >= 0:
        i += 1
    return i


def _hgvs_to_zbc(i: int):
    """Convert hgvs (1 based, missing zero)"""
    if i >= 1:
        i -= 1
    return i


class AlignmentMapper:
    """Maps hgvs location objects between genomic (g), non-coding (n) and
    cds (c) coordinates according to a CIGAR string.

    :param hdp: HGVS Data Provider Interface-compliant instance (see :class:`hgvs.dataproviders.interface.Interface`)
    :param str tx_ac: string representing transcript accession (e.g., NM_000551.2)
    :param str alt_ac: string representing the reference sequence accession (e.g., NC_000019.10)
    :param str alt_aln_method: string representing the alignment method; valid values depend on data source

    """

    __slots__ = (
        "tx_ac",
        "alt_ac",
        "alt_aln_method",
        "strand",
        "gc_offset",
        "cds_start_i",
        "cds_end_i",
        "tgt_len",
        "cigarmapper",
        "ref_pos",
        "tgt_pos",
        "cigar_op",
    )

    def __init__(self, hdp, tx_ac, alt_ac, alt_aln_method):
        self.tx_ac = tx_ac
        self.alt_ac = alt_ac
        self.alt_aln_method = alt_aln_method

        if self.alt_aln_method != "transcript":
            tx_info = hdp.get_tx_info(self.tx_ac, self.alt_ac, self.alt_aln_method)
            if tx_info is None:
                raise HGVSDataNotAvailableError(
                    "AlignmentMapper(tx_ac={self.tx_ac}, "
                    "alt_ac={self.alt_ac}, alt_aln_method={self.alt_aln_method}): "
                    "No transcript info".format(self=self)
                )

            tx_exons = hdp.get_tx_exons(self.tx_ac, self.alt_ac, self.alt_aln_method)
            if tx_exons is None:
                raise HGVSDataNotAvailableError(
                    "AlignmentMapper(tx_ac={self.tx_ac}, "
                    "alt_ac={self.alt_ac}, alt_aln_method={self.alt_aln_method}): "
                    "No transcript exons".format(self=self)
                )

            # hgvs-386: An assumption when building the cigar string
            # is that exons are adjacent. Assert that here.
            sorted_tx_exons = sorted(tx_exons, key=lambda e: e["ord"])
            for i in range(1, len(sorted_tx_exons)):
                if (
                    sorted_tx_exons[i - 1]["tx_end_i"]
                    != sorted_tx_exons[i]["tx_start_i"]
                ):
                    raise HGVSDataNotAvailableError(
                        "AlignmentMapper(tx_ac={self.tx_ac}, "
                        "alt_ac={self.alt_ac}, alt_aln_method={self.alt_aln_method}): "
                        "Exons {a} and {b} are not adjacent".format(
                            self=self, a=i, b=i + 1
                        )
                    )

            self.strand = tx_exons[0]["alt_strand"]
            self.gc_offset = tx_exons[0]["alt_start_i"]
            self.cds_start_i = tx_info["cds_start_i"]
            self.cds_end_i = tx_info["cds_end_i"]

            cigar = build_tx_cigar(tx_exons, self.strand)
            self.cigarmapper = CIGARMapper(cigar)
            self.tgt_len = self.cigarmapper.tgt_len

        else:
            # this covers the identity cases n <-> c
            tx_identity_info = hdp.get_tx_identity_info(self.tx_ac)
            if tx_identity_info is None:
                raise HGVSDataNotAvailableError(
                    "AlignmentMapper(tx_ac={self.tx_ac}, "
                    "alt_ac={self.alt_ac}, alt_aln_method={self.alt_aln_method}): "
                    "No transcript info".format(self=self)
                )
            self.cds_start_i = tx_identity_info["cds_start_i"]
            self.cds_end_i = tx_identity_info["cds_end_i"]
            self.tgt_len = sum(tx_identity_info["lengths"])
            self.cigarmapper = None

        assert not ((self.cds_start_i is None) ^ (self.cds_end_i is None)), (
            "CDS start and end must both be defined or neither defined"
        )

    def __str__(self):
        return (
            "{self.__class__.__name__}: {self.tx_ac} ~ {self.alt_ac} ~ {self.alt_aln_method}; "
            "{strand_pm} strand; offset={self.gc_offset}".format(
                self=self, strand_pm=strand_int_to_pm(self.strand)
            )
        )

    def _extract_genomic_position(
        self, interval_part, is_start: bool
    ) -> int | Interval:
        """Extract genomic position from an interval part (start or end).

        Args:
            interval_part: The start or end part of an interval
            is_start: True if this is the start part, False if end part

        Returns:
            int: The genomic position (zero-based relative to alignment)
        """
        if isinstance(interval_part, Interval):
            # For interval parts, use the appropriate boundary based on strand
            if self.strand == -1:
                # Reverse strand: use start for start, end for end
                if is_start:
                    return interval_part.start.base - 1 - self.gc_offset
                else:
                    return interval_part.end.base - 1 - self.gc_offset
            else:
                # Forward strand: use end for start, start for end
                if is_start:
                    return interval_part.end.base - 1 - self.gc_offset
                else:
                    return interval_part.start.base - 1 - self.gc_offset
        elif hasattr(interval_part, "start") and hasattr(interval_part, "end"):
            # Handle case where interval_part itself is an interval (like in _g_to_n_interval)
            # Extract start and end positions using _get_start_end
            start_pos, end_pos = self._get_start_end(interval_part)
            if is_start:
                return start_pos.base - 1 - self.gc_offset
            else:
                return end_pos.base - 1 - self.gc_offset

        else:
            # Simple position
            return interval_part.base - 1 - self.gc_offset

    def _get_start_end(self, var):
        """Get start and end positions from the variant.

        Args:
            var: A variant object with posedit.pos attribute

        Returns:
            tuple: (start_position, end_position) where positions can be SimplePosition or BaseOffsetPosition
        """
        from hgvs.utils.position import get_start_end

        return get_start_end(var)

    def _create_base_offset_interval(
        self,
        start_pos,
        start_offset,
        start_cigar,
        end_pos,
        end_offset,
        end_cigar,
        uncertain,
    ):
        """Create a BaseOffsetInterval from position data.

        Args:
            start_pos: Start position value
            start_offset: Start offset value
            start_cigar: Start cigar operation
            end_pos: End position value
            end_offset: End offset value
            end_cigar: End cigar operation
            uncertain: Whether the interval is uncertain

        Returns:
            BaseOffsetInterval: The created interval
        """

        if start_pos:
            start = hgvs.location.BaseOffsetPosition(
                base=_zbc_to_hgvs(start_pos),
                offset=start_offset,
                datum=Datum.SEQ_START,
                uncertain=start_cigar in "DI",
            )
        else:
            start = hgvs.location.BaseOffsetPosition(
                base=None,
                offset=0,
                datum=Datum.SEQ_START,
                uncertain=False,
            )

        if end_pos:
            end = hgvs.location.BaseOffsetPosition(
                base=_zbc_to_hgvs(end_pos),
                offset=end_offset,
                datum=Datum.SEQ_START,
                uncertain=end_cigar in "DI",
            )
        else:
            end = hgvs.location.BaseOffsetPosition(
                base=None,
                offset=0,
                datum=Datum.SEQ_START,
                uncertain=False,
            )

        return hgvs.location.BaseOffsetInterval(
            start=start,
            end=end,
            uncertain=uncertain,
        )

    def _g_to_n_interval(
        self, g_interval: Interval, strict_bounds: Optional[bool] = None
    ) -> BaseOffsetInterval:
        """Convert a genomic (g.) interval to a transcript cDNA (n.) interval.

        This method handles cases where both start and end of the input interval are themselves intervals.
        The returned interval is corrected based on the strand direction.

        Args:
            g_interval: A genomic interval where both start and end are intervals
            strict_bounds: Whether to enforce strict bounds checking

        Returns:
            A transcript cDNA interval with proper strand correction
        """
        # Get start and end positions from the genomic interval
        start_start, start_end = self._get_start_end(g_interval.start)
        end_start, end_end = self._get_start_end(g_interval.end)

        # Convert to zero-based coordinates relative to the alignment
        if not start_start.base:
            grs_start = None
        else:
            grs_start = start_start.base - 1 - self.gc_offset
        grs_end = start_end.base - 1 - self.gc_offset
        gre_start = end_start.base - 1 - self.gc_offset
        if not end_end.base:
            gre_end = None
        else:
            gre_end = end_end.base - 1 - self.gc_offset

        # Map genomic positions to transcript positions
        if grs_start:
            (
                n_start_start,
                n_start_start_offset,
                n_start_start_cigar,
            ) = self.cigarmapper.map_ref_to_tgt(
                pos=grs_start, end="start", strict_bounds=strict_bounds
            )
        else:
            n_start_start = None
            n_start_start_offset = 0
            n_start_start_cigar = None

        (
            n_start_end,
            n_start_end_offset,
            n_start_end_cigar,
        ) = self.cigarmapper.map_ref_to_tgt(
            pos=grs_end, end="end", strict_bounds=strict_bounds
        )
        (
            n_end_start,
            n_end_start_offset,
            n_end_start_cigar,
        ) = self.cigarmapper.map_ref_to_tgt(
            pos=gre_start, end="start", strict_bounds=strict_bounds
        )
        if gre_end:
            (
                n_end_end,
                n_end_end_offset,
                n_end_end_cigar,
            ) = self.cigarmapper.map_ref_to_tgt(
                pos=gre_end, end="end", strict_bounds=strict_bounds
            )
        else:
            n_end_end = None
            n_end_end_offset = 0
            n_end_end_cigar = None

        # For reverse strand, transform coordinates
        if self.strand == -1:
            # Store original values
            orig_n_start_start = n_start_start
            orig_n_start_end = n_start_end
            orig_n_end_start = n_end_start
            orig_n_end_end = n_end_end

            # Transform positions
            n_start_start = self.tgt_len - 1 - orig_n_start_end
            if orig_n_start_start:
                n_start_end = self.tgt_len - 1 - orig_n_start_start
            else:
                n_start_end = None
            if orig_n_end_end:
                n_end_start = self.tgt_len - 1 - orig_n_end_end
            else:
                n_end_start = None
            n_end_end = self.tgt_len - 1 - orig_n_end_start

            # Transform offsets
            n_start_start_offset, n_start_end_offset = (
                -n_start_end_offset,
                -n_start_start_offset,
            )
            n_end_start_offset, n_end_end_offset = (
                -n_end_end_offset,
                -n_end_start_offset,
            )

            # Swap cigar operations
            n_start_start_cigar, n_start_end_cigar = (
                n_start_end_cigar,
                n_start_start_cigar,
            )
            n_end_start_cigar, n_end_end_cigar = n_end_end_cigar, n_end_start_cigar

        # Create the start interval
        start_interval = self._create_base_offset_interval(
            n_start_start,
            n_start_start_offset,
            n_start_start_cigar,
            n_start_end,
            n_start_end_offset,
            n_start_end_cigar,
            g_interval.start.uncertain,
        )

        # Create the end interval
        end_interval = self._create_base_offset_interval(
            n_end_start,
            n_end_start_offset,
            n_end_start_cigar,
            n_end_end,
            n_end_end_offset,
            n_end_end_cigar,
            g_interval.end.uncertain,
        )

        # For reverse strand, ensure start is less than or equal to end
        if self.strand == -1:
            if (
                start_interval.end.base
                and end_interval.start.base
                and start_interval.end.base > end_interval.start.base
            ):
                start_interval, end_interval = end_interval, start_interval
            elif (
                start_interval.start.base
                and end_interval.end.base
                and start_interval.start.base > end_interval.end.base
            ):
                start_interval, end_interval = end_interval, start_interval

        # Return the final interval
        return hgvs.location.Interval(
            start=start_interval,
            end=end_interval,
            uncertain=g_interval.uncertain,
        )

    def g_to_n(
        self, g_interval: Interval, strict_bounds: Optional[bool] = None
    ) -> BaseOffsetInterval:
        """convert a genomic (g.) interval to a transcript cDNA (n.) interval"""

        if strict_bounds is None:
            strict_bounds = global_config.mapping.strict_bounds

        if (
            isinstance(g_interval.start, Interval)
            # and g_interval.start.start.base
            and isinstance(g_interval.end, Interval)
            # and g_interval.end.end.base
        ):
            return self._g_to_n_interval(g_interval, strict_bounds)

        # in case of uncertain ranges, we fall back to the inner (more confident) interval
        grs = self._extract_genomic_position(g_interval.start, True)
        gre = self._extract_genomic_position(g_interval.end, False)

        forward_rna_start, forward_rna_start_offset, forward_rna_start_cigar = (
            self.cigarmapper.map_ref_to_tgt(
                pos=grs, end="start", strict_bounds=strict_bounds
            )
        )
        forward_rna_end, forward_rna_end_offset, forward_rna_end_cigar = (
            self.cigarmapper.map_ref_to_tgt(
                pos=gre, end="end", strict_bounds=strict_bounds
            )
        )

        if self.strand == -1:
            # Store original values before transformation
            orig_start = forward_rna_start
            orig_end = forward_rna_end

            # Transform coordinates for reverse strand
            forward_rna_start = self.tgt_len - 1 - orig_end
            forward_rna_end = self.tgt_len - 1 - orig_start

            orig_offset_start = forward_rna_start_offset
            orig_offset_end = forward_rna_end_offset
            forward_rna_start_offset = -orig_offset_end
            forward_rna_end_offset = -orig_offset_start

        start = hgvs.location.BaseOffsetPosition(
            base=_zbc_to_hgvs(forward_rna_start),
            offset=forward_rna_start_offset,
            datum=Datum.SEQ_START,
            uncertain=g_interval.start.uncertain,
        )
        end = hgvs.location.BaseOffsetPosition(
            base=_zbc_to_hgvs(forward_rna_end),
            offset=forward_rna_end_offset,
            datum=Datum.SEQ_START,
            uncertain=g_interval.end.uncertain,
        )

        # The returned interval would be uncertain when locating at alignment gaps
        # of if the initial interval was uncertain
        final_interval = hgvs.location.BaseOffsetInterval(
            start=start,
            end=end,
            uncertain=forward_rna_start_cigar in "DI" or forward_rna_end_cigar in "DI",
        )
        return final_interval

    def n_to_g(
        self, n_interval: Interval, strict_bounds: Optional[bool] = None
    ) -> Interval:
        """Convert a transcript (n.) interval to a genomic (g.) interval.

        Args:
            n_interval: A transcript interval
            strict_bounds: Whether to enforce strict bounds checking

        Returns:
            A genomic interval
        """
        if strict_bounds is None:
            strict_bounds = global_config.mapping.strict_bounds

        if isinstance(n_interval.start, Interval) and isinstance(
            n_interval.end, Interval
        ):
            return self._n_to_g_interval(n_interval, strict_bounds)

        # Get start and end positions
        s, e = self._get_start_end(n_interval)

        # Convert to zero-based coordinates
        frs = _hgvs_to_zbc(s.base)
        start_offset = s.offset
        fre = _hgvs_to_zbc(e.base)
        end_offset = e.offset

        # For reverse strand, transform coordinates before mapping
        if self.strand == -1:
            # Store original values
            orig_frs = frs
            orig_fre = fre
            orig_start_offset = start_offset
            orig_end_offset = end_offset

            # Transform coordinates
            frs = self.tgt_len - 1 - orig_fre
            fre = self.tgt_len - 1 - orig_frs
            start_offset = -orig_end_offset
            end_offset = -orig_start_offset

        # Map transcript positions to genomic positions
        grs, _, grs_cigar = self.cigarmapper.map_tgt_to_ref(
            pos=frs, end="start", strict_bounds=strict_bounds
        )
        gre, _, gre_cigar = self.cigarmapper.map_tgt_to_ref(
            pos=fre, end="end", strict_bounds=strict_bounds
        )

        # Add offset to get final genomic positions
        grs = grs + self.gc_offset + 1
        gre = gre + self.gc_offset + 1
        gs = grs + start_offset
        ge = gre + end_offset

        # Handle uncertain positions
        if isinstance(n_interval.start, Interval):
            start = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(gs, uncertain=False),
                end=hgvs.location.SimplePosition(uncertain=False),
                uncertain=True,
            )
        else:
            start = hgvs.location.SimplePosition(gs, uncertain=False)

        if isinstance(n_interval.end, Interval):
            end = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(ge, uncertain=False),
                end=hgvs.location.SimplePosition(uncertain=False),
                uncertain=True,
            )
        else:
            end = hgvs.location.SimplePosition(ge, uncertain=False)

        if isinstance(start, Interval) and isinstance(end, Interval):
            return hgvs.location.Interval(
                start=start,
                end=end,
                uncertain=grs_cigar in "DI" or gre_cigar in "DI",
            )

        # For reverse strand, ensure start is less than or equal to end
        if self.strand == -1:
            if start.base > end.base:
                start, end = end, start

        # The returned interval would be uncertain when locating at alignment gaps
        return hgvs.location.Interval(
            start=start,
            end=end,
            uncertain=grs_cigar in "DI" or gre_cigar in "DI",
        )

    def n_to_c(
        self,
        n_interval: Interval,
        strict_bounds: bool | None = None,
    ):
        """convert a transcript cDNA (n.) interval to a transcript CDS (c.) interval"""

        if strict_bounds is None:
            strict_bounds = global_config.mapping.strict_bounds

        if (
            self.cds_start_i is None
        ):  # cds_start_i defined iff cds_end_i defined; see assertion above
            raise HGVSUsageError(
                "CDS is undefined for {self.tx_ac}; cannot map to c. coordinate (non-coding transcript?)".format(
                    self=self
                )
            )

        if isinstance(n_interval, BaseOffsetInterval):
            start = n_interval.start.base
            end = n_interval.end.base
        elif isinstance(n_interval, Interval):
            start = n_interval.start.start.base
            end = n_interval.end.end.base

        if not start:
            start = 1

        if strict_bounds and (start <= 0 or (end and end > self.tgt_len)):
            raise HGVSInvalidIntervalError(
                "The given coordinate is outside the bounds of the reference sequence."
            )

        def pos_n_to_c(pos) -> hgvs.location.BaseOffsetPosition:
            if pos.base <= self.cds_start_i:
                c = pos.base - self.cds_start_i - (1 if pos.base > 0 else 0)
                c_datum = Datum.CDS_START
            elif pos.base > self.cds_start_i and pos.base <= self.cds_end_i:
                c = pos.base - self.cds_start_i
                c_datum = Datum.CDS_START
            else:
                c = pos.base - self.cds_end_i
                c_datum = Datum.CDS_END
            return hgvs.location.BaseOffsetPosition(
                base=c, offset=pos.offset, datum=c_datum, uncertain=pos.uncertain
            )

        def _interval_n_to_c(interval: Interval) -> hgvs.location.BaseOffsetInterval:
            istart = interval.start
            if istart and istart.base:
                c_start = pos_n_to_c(istart)
            else:
                c_start = hgvs.location.BaseOffsetPosition(
                    base=None, datum=Datum.CDS_START
                )
                if not interval.end or not interval.end.base:
                    raise HGVSInvalidIntervalError(
                        "n_to_c: interval start is None and end is None"
                    )

            iend = interval.end
            if iend and iend.base:
                c_end = pos_n_to_c(iend)
            else:
                c_end = hgvs.location.BaseOffsetPosition(base=None, datum=Datum.CDS_END)

            return hgvs.location.BaseOffsetInterval(
                start=c_start,
                end=c_end,
                uncertain=interval.uncertain,
            )

        baseoffset = True
        if isinstance(n_interval.start, Interval):
            c_start = _interval_n_to_c(n_interval.start)
            baseoffset = False
        else:
            c_start = pos_n_to_c(n_interval.start)

        if isinstance(n_interval.end, Interval):
            c_end = _interval_n_to_c(n_interval.end)
            baseoffset = False
        else:
            # pos
            c_end = pos_n_to_c(n_interval.end)

        if baseoffset:
            c_interval = hgvs.location.BaseOffsetInterval(
                start=c_start,
                end=c_end,
                uncertain=n_interval.uncertain,
            )
        else:
            c_interval = hgvs.location.Interval(
                start=c_start,
                end=c_end,
                uncertain=n_interval.uncertain,
            )
        return c_interval

    def c_to_n(self, c_interval, strict_bounds=None):
        """convert a transcript CDS (c.) interval to a transcript cDNA (n.) interval"""

        if strict_bounds is None:
            strict_bounds = global_config.mapping.strict_bounds

        if self.cds_start_i is None:
            raise HGVSUsageError(
                "CDS is undefined for {self.tx_ac}; this accession appears to be for a non-coding transcript".format(
                    self=self
                )
            )

        def pos_c_to_n(pos):
            if not pos or not pos.base:
                return hgvs.location.BaseOffsetPosition(
                    base=None, datum=Datum.SEQ_START
                )

            if pos.datum == Datum.CDS_START:
                n = pos.base + self.cds_start_i
                if pos.base < 0:  # correct for lack of c.0 coordinate
                    n += 1
            elif pos.datum == Datum.CDS_END:
                n = pos.base + self.cds_end_i
            if n <= 0:  # correct for lack of n.0 coordinate
                n -= 1
            if n <= 0 or n > self.tgt_len:
                if strict_bounds:
                    raise HGVSInvalidIntervalError(
                        f"c.{pos} coordinate is out of bounds"
                    )

            return hgvs.location.BaseOffsetPosition(
                base=n,
                offset=pos.offset,
                datum=Datum.SEQ_START,
                uncertain=pos.uncertain,
            )

        def _interval_c_to_n(interval: Interval) -> hgvs.location.BaseOffsetInterval:
            istart = interval.start
            n_start = pos_c_to_n(istart)

            iend = interval.end
            n_end = pos_c_to_n(iend)

            return hgvs.location.BaseOffsetInterval(
                start=n_start,
                end=n_end,
                uncertain=interval.uncertain,
            )

        baseoffset = True
        if isinstance(c_interval.start, Interval):
            n_start = _interval_c_to_n(c_interval.start)
            baseoffset = False
        else:
            n_start = pos_c_to_n(c_interval.start)

        if isinstance(c_interval.end, Interval):
            n_end = _interval_c_to_n(c_interval.end)
            baseoffset = False
        else:
            n_end = pos_c_to_n(c_interval.end)

        if baseoffset:
            n_interval = hgvs.location.BaseOffsetInterval(
                start=n_start,
                end=n_end,
                uncertain=c_interval.uncertain,
            )
        else:
            n_interval = hgvs.location.Interval(
                start=n_start,
                end=n_end,
                uncertain=c_interval.uncertain,
            )
        return n_interval

    def g_to_c(self, g_interval, strict_bounds: bool | None = None):
        """convert a genomic (g.) interval to a transcript CDS (c.) interval"""
        var_n = self.g_to_n(g_interval)
        return self.n_to_c(var_n, strict_bounds=strict_bounds)

    def c_to_g(self, c_interval, strict_bounds=None):
        """convert a transcript CDS (c.) interval to a genomic (g.) interval"""

        var_n = self.c_to_n(c_interval)
        return self.n_to_g(var_n, strict_bounds=strict_bounds)

    @property
    def is_coding_transcript(self):
        if (self.cds_start_i is not None) ^ (self.cds_end_i is not None):
            raise HGVSError(
                "{self.tx_ac}: CDS start_i and end_i"
                " must be both defined or both undefined".format(self=self)
            )
        return self.cds_start_i is not None

    def g_interval_is_inbounds(self, ival):
        grs = ival.start.base - 1 - self.gc_offset
        gre = ival.end.base - 1 - self.gc_offset
        return grs >= 0 and gre <= self.cigarmapper.ref_len

    def fix_offset(self, pos):
        if isinstance(pos, int):
            return pos
        if not pos.base:
            return None

        return _hgvs_to_zbc(pos.base)

    def _n_to_g_interval(
        self, n_interval: Interval, strict_bounds: Optional[bool] = None
    ) -> Interval:
        """Convert transcript (n.) intervals to a genomic (g.) interval.

        This method handles cases where both start and end of the input interval are themselves intervals.
        The returned interval is corrected based on the strand direction.

        Args:
            n_interval: A transcript interval where both start and end are intervals
            strict_bounds: Whether to enforce strict bounds checking

        Returns:
            A genomic interval with proper strand correction
        """
        # Get start and end positions from the intervals
        start_start, start_end = self._get_start_end(n_interval.start)
        end_start, end_end = self._get_start_end(n_interval.end)

        # Convert to zero-based coordinates
        frs_start = self.fix_offset(start_start)
        frs_end = self.fix_offset(start_end)
        fre_start = self.fix_offset(end_start)
        fre_end = self.fix_offset(end_end)

        ss_offset = start_start.offset
        se_offset = start_end.offset
        es_offset = end_start.offset
        ee_offset = end_end.offset
        # For reverse strand, transform coordinates before mapping
        if self.strand == -1:
            # Store original values
            orig_frs_start = frs_start
            orig_frs_end = frs_end
            orig_fre_start = fre_start
            orig_fre_end = fre_end
            orig_ss_offset = ss_offset
            orig_se_offset = se_offset
            orig_es_offset = es_offset
            orig_ee_offset = ee_offset

            # Transform coordinates
            frs_start = self.tgt_len - 1 - orig_frs_end
            if orig_frs_start:
                frs_end = self.tgt_len - 1 - orig_frs_start
            else:
                frs_end = None
            fre_start = self.tgt_len - 1 - orig_fre_start
            if orig_fre_end:
                fre_end = self.tgt_len - 1 - orig_fre_end
            else:
                fre_end = None
            ss_offset = orig_ss_offset
            se_offset = orig_se_offset
            es_offset = orig_es_offset
            ee_offset = orig_ee_offset

        left_boundary_out_of_bounds = False
        right_boundary_out_of_bounds = False

        # Map transcript positions to genomic positions
        if frs_start is not None:
            grs_start, _, grs_start_cigar = self.cigarmapper.map_tgt_to_ref(
                pos=frs_start, end="start", strict_bounds=strict_bounds
            )
        else:
            grs_start = None
            grs_start_cigar = None

        if frs_end is not None:
            grs_end, _, grs_end_cigar = self.cigarmapper.map_tgt_to_ref(
                pos=frs_end, end="end", strict_bounds=strict_bounds
            )
        else:
            grs_end = None
            grs_end_cigar = None

        # is gre_start ever None?
        gre_start, _, gre_start_cigar = self.cigarmapper.map_tgt_to_ref(
            pos=fre_start, end="start", strict_bounds=strict_bounds
        )
        if fre_end is not None:
            gre_end, _, gre_end_cigar = self.cigarmapper.map_tgt_to_ref(
                pos=fre_end, end="end", strict_bounds=strict_bounds
            )
        else:
            gre_end = None
            gre_end_cigar = None

        # Add offset to get final genomic positions
        if self.strand == -1:
            if gre_start is not None:
                gre_orig_start = gre_start + self.gc_offset + 1 - es_offset
            else:
                gre_orig_start = -1

            if gre_end is not None:
                gre_orig_end = gre_end + self.gc_offset + 1 - ee_offset
            else:
                gre_orig_end = -1
            if grs_start is not None:
                grs_orig_start = grs_start + self.gc_offset + 1 - se_offset
            else:
                grs_orig_start = math.inf
            if grs_end is not None:
                grs_orig_end = grs_end + self.gc_offset + 1 - ss_offset
            else:
                grs_orig_end = math.inf + 1

            # Sort the four values from smallest to largest
            all_positions = [gre_orig_start, gre_orig_end, grs_orig_start, grs_orig_end]
            all_positions.sort()

            # Convert sentinel values back to None
            all_positions = [
                None if x == -1 or x >= math.inf else x for x in all_positions
            ]

            grs_start, grs_end, gre_start, gre_end = all_positions

        else:
            if grs_start:
                grs_start = grs_start + self.gc_offset + 1 + ss_offset
            if grs_end:
                grs_end = grs_end + self.gc_offset + 1 + se_offset
            if gre_start:
                gre_start = gre_start + self.gc_offset + 1 + es_offset
            if gre_end:
                gre_end = gre_end + self.gc_offset + 1 + ee_offset

        # For reverse strand, ensure start is less than or equal to end

        if left_boundary_out_of_bounds:
            lstart = hgvs.location.SimplePosition(
                base=None, uncertain=start_start.uncertain
            )
        else:
            lstart = hgvs.location.SimplePosition(
                base=grs_start, uncertain=start_start.uncertain
            )

        lend = hgvs.location.SimplePosition(grs_end, uncertain=start_end.uncertain)

        rstart = hgvs.location.SimplePosition(gre_start, uncertain=end_start.uncertain)
        if right_boundary_out_of_bounds:
            rend = hgvs.location.SimplePosition(
                base=None, uncertain=end_start.uncertain
            )
        else:
            rend = hgvs.location.SimplePosition(
                base=gre_end, uncertain=end_start.uncertain
            )

        # Create the start interval
        g_start = hgvs.location.Interval(
            start=lstart,
            end=lend,
            uncertain=n_interval.start.uncertain
            or grs_start_cigar in "DI"
            or grs_end_cigar in "DI",
        )

        # Create the end interval
        g_end = hgvs.location.Interval(
            start=rstart,
            end=rend,
            uncertain=n_interval.end.uncertain
            or gre_start_cigar in "DI"
            or gre_end_cigar in "DI",
        )

        if self.strand == -1:
            if (
                g_start.end.base
                and g_end.start.base
                and g_start.end.base > g_end.start.base
            ):
                g_start, g_end = g_end, g_start

        # Return the final interval
        return hgvs.location.Interval(
            start=g_start,
            end=g_end,
            uncertain=n_interval.uncertain,
        )


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
