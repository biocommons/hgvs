"""
Utility class that projects variants from one transcript to another via a
common reference sequence.
"""

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

    #TODO: project_variant_{forward,backward}
