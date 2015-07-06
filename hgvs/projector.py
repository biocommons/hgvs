# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
"""
Utility class that projects variants from one transcript to another via a
common reference sequence.
"""

import copy

import hgvs.transcriptmapper


class Projector(object):
    """
    The Projector class implements liftover between two transcripts via a
    common reference sequence.

    :param hdp: HGVS Data Provider Interface-compliant instance (see :class:`hgvs.dataproviders.interface.Interface`)
    :param ref: string representing the common reference assembly (e.g., GRCh37.p10)
    :param src_ac: string representing the source transcript accession (e.g., NM_000551.2)
    :param dst_ac: string representing the destination transcript accession (e.g., NM_000551.3)
    :param src_alt_aln_method: string representing the source transcript alignment method
    :param dst_alt_aln_method: string representing the destination transcript alignment method

    This class assumes (and verifies) that the transcripts are on the same
    strand. This assumption obviates some work in flipping sequence
    variants twice unnecessarily.
    """

    def __init__(self, hdp, alt_ac, src_ac, dst_ac, src_alt_aln_method='splign', dst_alt_aln_method='splign'):
        self.hdp = hdp
        self.alt_ac = alt_ac
        self.src_tm = hgvs.transcriptmapper.TranscriptMapper(hdp, src_ac, alt_ac, src_alt_aln_method)
        self.dst_tm = hgvs.transcriptmapper.TranscriptMapper(hdp, dst_ac, alt_ac, dst_alt_aln_method)

    def project_interval_forward(self, c_interval):
        """
        project c_interval on the source transcript to the
        destination transcript

        :param c_interval: an :class:`hgvs.interval.Interval` object on the source transcript
        :returns: c_interval: an :class:`hgvs.interval.Interval` object on the destination transcript
        """
        return self.dst_tm.g_to_c(self.src_tm.c_to_g(c_interval))

    def project_interval_backward(self, c_interval):
        """
        project c_interval on the destination transcript to the
        source transcript

        :param c_interval: an :class:`hgvs.interval.Interval` object on the destination transcript
        :returns: c_interval: an :class:`hgvs.interval.Interval` object on the source transcript
        """
        return self.src_tm.g_to_c(self.dst_tm.c_to_g(c_interval))

    def project_variant_forward(self, c_variant):
        """
        project c_variant on the source transcript onto the destination transcript

        :param c_variant: an :class:`hgvs.variant.SequenceVariant` object on the source transcript
        :returns: c_variant: an :class:`hgvs.variant.SequenceVariant` object on the destination transcript
        """
        if c_variant.ac != self.src_tm.tx_ac:
            raise RuntimeError('variant accession does not match that used to initialize ' + __name__)
        new_c_variant = copy.deepcopy(c_variant)
        new_c_variant.ac = self.dst_tm.tx_ac
        new_c_variant.posedit.pos = self.project_interval_forward(c_variant.posedit.pos)
        return new_c_variant

    def project_variant_backward(self, c_variant):
        """
        project c_variant on the source transcript onto the destination transcript

        :param c_variant: an :class:`hgvs.variant.SequenceVariant` object on the source transcript
        :returns: c_variant: an :class:`hgvs.variant.SequenceVariant` object on the destination transcript
        """
        if c_variant.ac != self.dst_tm.tx_ac:
            raise RuntimeError('variant accession does not match that used to initialize ' + __name__)
        new_c_variant = copy.deepcopy(c_variant)
        new_c_variant.ac = self.src_tm.tx_ac
        new_c_variant.posedit.pos = self.project_interval_backward(c_variant.posedit.pos)
        return new_c_variant
