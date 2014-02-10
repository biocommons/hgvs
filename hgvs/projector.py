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

    :param bdi: Bioinformatics Data Interface-compliant instance (see :class:`bdi.interface0.Interface0`)
    :param ref: string representing the common reference assembly (e.g., GRCh37.p10)
    :param src_ac: string representing the source transcript accession (e.g., NM_000551.2)
    :param dst_ac: string representing the destination transcript accession (e.g., NM_000551.3)

    This class assumes (and verifies) that the transcripts are on the same
    strand. This assumption obviates some work in flipping sequence
    variants twice unnecessarily.
    """
    def __init__(self,bdi,ref,src_ac,dst_ac):
        self.bdi = bdi
        self.ref = ref
        self.src_tm = hgvs.transcriptmapper.TranscriptMapper(bdi,src_ac,ref)
        self.dst_tm = hgvs.transcriptmapper.TranscriptMapper(bdi,dst_ac,ref)

    def project_interval_forward(self,c_interval):
        """
        project c_interval on the source transcript to the
        destination transcript

        :param c_interval: an :class:`hgvs.interval.Interval` object on the source transcript
        :returns: c_interval: an :class:`hgvs.interval.Interval` object on the destination transcript
        """
        return self.dst_tm.hgvsg_to_hgvsc( self.src_tm.hgvsc_to_hgvsg( c_interval ) )
        
    def project_interval_backward(self,c_interval):
        """
        project c_interval on the destination transcript to the
        source transcript

        :param c_interval: an :class:`hgvs.interval.Interval` object on the destination transcript
        :returns: c_interval: an :class:`hgvs.interval.Interval` object on the source transcript
        """
        return self.src_tm.hgvsg_to_hgvsc( self.dst_tm.hgvsc_to_hgvsg( c_interval ) )


    def project_variant_forward(self,c_variant):
        """
        project c_variant on the source transcript onto the destination transcript

        :param c_variant: an :class:`hgvs.variant.SequenceVariant` object on the source transcript
        :returns: c_variant: an :class:`hgvs.variant.SequenceVariant` object on the destination transcript
        """
        if c_variant.ac != self.src_tm.ac:
            raise RuntimeError('variant accession does not match that used to initialize '+__name__)
        new_c_variant = copy.deepcopy( c_variant )
        new_c_variant.ac = self.dst_tm.ac
        new_c_variant.posedit.pos = self.project_interval_forward( c_variant.posedit.pos )
        return new_c_variant

    def project_variant_backward(self,c_variant):
        """
        project c_variant on the source transcript onto the destination transcript

        :param c_variant: an :class:`hgvs.variant.SequenceVariant` object on the source transcript
        :returns: c_variant: an :class:`hgvs.variant.SequenceVariant` object on the destination transcript
        """
        if c_variant.ac != self.dst_tm.ac:
            raise RuntimeError('variant accession does not match that used to initialize '+__name__)
        new_c_variant = copy.deepcopy( c_variant )
        new_c_variant.ac = self.src_tm.ac
        new_c_variant.posedit.pos = self.project_interval_backward( c_variant.posedit.pos )
        return new_c_variant


