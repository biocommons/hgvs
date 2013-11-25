import hgvs.edti.interface

class UTA0(hgvs.edti.interface.Interface):
    def __init__(self,uta_conn):
        self.uta_conn = uta_conn

    def fetch_gene_info(self,gene):
        pass

    def fetch_gene_transcripts(self,gene):
        return self.uta_conn.get_tx_for_gene(gene)

    def fetch_transcript_info(self,ac):
        return self.uta_conn.get_tx_info(ac)

    def fetch_transcript_exons(self,ac,assy):
        return self.uta_conn.get_tx_exons(ac,assy)
