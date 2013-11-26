import recordtype


class TranscriptInfo( recordtype.recordtype('TranscriptInfo', [
        'ac', 'gene', 'chrom', 'strand', 'cds_start_i', 'cds_end_i',
        'gene_description', 'gene_summary' ])):
    pass

class TranscriptExon( recordtype.recordtype('TranscriptExon', [
    'ac', 'assy', 'ord', 'name', 'start_i', 'end_i', 'g_start_i', 'g_end_i', 'cigar', ])):
    pass




    def get_tx_info(self,ac):
        """return transcript info for supplied accession (ac), or None if not found

        :param ac: transcript accession with version (e.g., 'NM_000051.3')
        :type ac: str
        """

        self.cur.execute(self.tx_info_sql,{'ac': ac})
        assert self.cur.rowcount <= 1, 'get_tx_info({ac}) unexpectedly returned {c} rows'.format(
            ac=ac, c=self.cur.rowcount)
        return self.cur.fetchone()        # None if no match

    def get_tx_exons(self,ac,ref):    
        """
        return transcript info for supplied accession (ac), or None if not found
        
        :param ac: transcript accession with version (e.g., 'NM_000051.3')
        :type ac: str
        :param ref: reference genome ('GRCh37.p10' is the only valid value at this time)
        :type ref: str
        
        # tx_exons = db.get_tx_exons('NM_000051.3','GRCh37.p10')
        # len(tx_exons)
        63
        
        tx_exons have the following attributes::
        
          {'ac': 'NM_000051.3',               # transcript accession
          'ref': 'GRCh37.p10',
          'g_start_i': 108093558,             # genomic start coordinate
          'g_end_i': 108093913,               # genome end coordinate
          'name': '1',
          'ord': 1,
          't_start_i': 0
          't_end_i': 355,
          'g_cigar': '355M',                  # CIGAR string, relative to genome
          'g_seq_a': None,
          't_seq_a': None,
          }
        
        .. note:: chromosome and strand are in the tx_info record
        
        For example:
        
        # tx_exons[0]['ac']
        'NM_000051.3'
        
        """
        self.cur.execute(self.tx_exons_sql,{'ac': ac, 'ref': ref})
        return self.cur.fetchall()        # [] if no match

    def get_tx_for_gene(self,gene):
        """
        return transcript info records for supplied gene, in order of decreasing length

        :param gene: HGNC gene name
        :type gene: str
        """
        self.cur.execute(self.tx_for_gene_sql,{'gene': gene})
        return self.cur.fetchall()
