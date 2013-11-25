from hgvs.edti.interface import Interface

class UTA0(Interface):
    def __init__(self,uta_conn):
        self.uta_conn = uta_conn

    def fetch_gene_info(self,gene):
        # Cheat: UTA doesn't have a fetch_gene_info method yet, so just
        # fetch transcripts and filter the result.
        gts = self.fetch_gene_transcripts(gene)
        if len(gts) == 0:
            return None
        gt = gts[0]
        return { k: gt[k] for k in ['gene','descr','summary','chr','strand'] }

    def fetch_gene_transcripts(self,gene):
        return self.uta_conn.get_tx_for_gene(gene)

    def fetch_transcript_info(self,ac):
        return self.uta_conn.get_tx_info(ac)

    def fetch_transcript_exons(self,ac,assy):
        return self.uta_conn.get_tx_exons(ac,assy)
