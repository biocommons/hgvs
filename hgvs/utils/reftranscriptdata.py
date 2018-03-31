from Bio.Seq import Seq

from hgvs.exceptions import HGVSDataNotAvailableError


class RefTranscriptData(object):
    def __init__(self, hdp, tx_ac, pro_ac):
        """helper for generating RefTranscriptData from for c_to_p"""
        tx_info = hdp.get_tx_identity_info(tx_ac)
        tx_seq = hdp.get_seq(tx_ac)

        if tx_info is None or tx_seq is None:
            raise HGVSDataNotAvailableError("Missing transcript data for accession: {}".format(tx_ac))

        # use 1-based hgvs coords
        cds_start = tx_info["cds_start_i"] + 1
        cds_stop = tx_info["cds_end_i"]

        # padding list so biopython won't complain during the conversion
        tx_seq_to_translate = tx_seq[cds_start - 1:cds_stop]
        if len(tx_seq_to_translate) % 3 != 0:
            "".join(list(tx_seq_to_translate).extend(["N"] * ((3 - len(tx_seq_to_translate) % 3) % 3)))

        tx_seq_cds = Seq(tx_seq_to_translate)
        protein_seq = str(tx_seq_cds.translate())

        if pro_ac is None:
            # get_acs... will always return at least the MD5_ accession
            pro_ac = (hdp.get_pro_ac_for_tx_ac(tx_ac) or hdp.get_acs_for_protein_seq(protein_seq)[0])

        self.transcript_sequence = tx_seq
        self.aa_sequence = protein_seq
        self.cds_start = cds_start
        self.cds_stop = cds_stop
        self.protein_accession = pro_ac
