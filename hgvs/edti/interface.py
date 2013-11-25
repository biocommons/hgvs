import abc

class Interface(object):
    __metaclass__ = abc.ABCMeta

    """Pure virtural class for the HGVS External Data and Tools Interface.
    Every EDTI implementation should be subclass (possibly indirect) of this class.
    """

    @abc.abstractmethod
    def fetch_gene_info(self,ac):
        pass

    @abc.abstractmethod
    def fetch_gene_transcripts(self,ac):
        pass

    @abc.abstractmethod
    def fetch_transcript_exons(self,ac):
        pass

    @abc.abstractmethod
    def fetch_transcript_info(self,ac):
        pass

#    @abc.abstractmethod
#    def fetch_chromosome_slice(self,chr,gs_i,ge_i):
#        pass
#
#    @abc.abstractmethod
#    def fetch_protein_ac_for_sequence(self,s):
#        pass
#
#    @abc.abstractmethod
#    def translate_mrna(self,s):
#        pass
