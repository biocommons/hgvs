import logging

from bioutils.sequences import TranslationTable, translate_cds

from hgvs.exceptions import HGVSDataNotAvailableError

# Constants for mitochondrial transcript handling
_REMAINDER_TWO = 2

_logger = logging.getLogger(__name__)


class RefTranscriptData:
    def __init__(self, hdp, tx_ac, pro_ac, translation_table=TranslationTable.standard):
        """helper for generating RefTranscriptData from for c_to_p"""
        tx_info = hdp.get_tx_identity_info(tx_ac)
        tx_seq = hdp.get_seq(tx_ac)

        if tx_info is None or tx_seq is None:
            msg = f"Missing transcript data for accession: {tx_ac}"
            raise HGVSDataNotAvailableError(msg)

        # use 1-based hgvs coords
        cds_start = tx_info["cds_start_i"] + 1
        cds_stop = tx_info["cds_end_i"]

        tx_seq_to_translate = tx_seq[cds_start - 1 : cds_stop]
        # Complete a partial stop codon (T or TA at the CDS end) with poly(A) bases.
        # This is the norm for mitochondrial transcripts, but RefSeq also annotates
        # nuclear transcripts whose stop codon is completed by polyadenylation
        # (e.g. NM_001424184.1, CDS length 659); ported from openvar/vv_hgvs.
        poly_a = self._poly_a_completion(tx_seq_to_translate)
        if poly_a:
            if translation_table != TranslationTable.vertebrate_mitochondrial:
                _logger.info(
                    "Completing partial stop codon of %s with poly(A) tail (%s)", tx_ac, poly_a
                )
            tx_seq = tx_seq[:cds_stop] + poly_a + tx_seq[cds_stop:]
            tx_seq_to_translate += poly_a
        # coding sequences that are not divisable by 3 are not yet supported
        if len(tx_seq_to_translate) % 3 != 0:
            msg = f"Transcript {tx_ac} is not supported because its sequence length of {len(tx_seq_to_translate)} is not divisible by 3."
            raise NotImplementedError(msg)

        protein_seq = translate_cds(tx_seq_to_translate, translation_table=translation_table)

        if pro_ac is None:
            # get_acs... will always return at least the MD5_ accession
            # TODO: drop get_acs_for_protein_seq; use known mapping or digest (wo/pro ac inference)
            pro_ac = hdp.get_pro_ac_for_tx_ac(tx_ac) or hdp.get_acs_for_protein_seq(protein_seq)[0]

        # Auto-detect selenocysteine: check if reference protein sequence contains 'U'
        # If it does, re-translate using the selenocysteine translation table
        actual_translation_table = translation_table
        try:
            ref_protein_seq = hdp.get_seq(pro_ac) if pro_ac else None
            if ref_protein_seq and "U" in ref_protein_seq:
                # Re-translate with selenocysteine table
                actual_translation_table = TranslationTable.selenocysteine
                protein_seq = translate_cds(
                    tx_seq_to_translate, translation_table=actual_translation_table
                )
        except Exception as e:
            # If we can't get reference sequence, use what we have
            # This is expected for some accessions that may not be available
            _logger.debug("Could not fetch reference protein sequence for %s: %s", pro_ac, e)

        exon_start_positions = [-tx_info["cds_start_i"]]
        exon_end_positions = [exon_start_positions[0] + tx_info["lengths"][0]]
        for exon_length in tx_info["lengths"][1:]:
            exon_start_positions.append(exon_end_positions[-1] + 1)
            exon_end_positions.append(exon_end_positions[-1] + exon_length)

        self.transcript_sequence = tx_seq
        self.aa_sequence = protein_seq
        self.cds_start = cds_start
        self.cds_stop = cds_stop
        self.protein_accession = pro_ac
        self.translation_table = actual_translation_table
        self.exon_start_positions = exon_start_positions
        self.exon_end_positions = exon_end_positions

    @staticmethod
    def _poly_a_completion(cds_seq):
        """Return the poly(A) bases that complete a partial stop codon at the end of
        cds_seq ("AA" for a trailing T, "A" for a trailing TA), or None if the CDS
        is already a whole number of codons or cannot be completed to a stop."""
        remainder = len(cds_seq) % 3
        if remainder == 1 and cds_seq.endswith("T"):
            return "AA"
        if remainder == _REMAINDER_TWO and cds_seq.endswith("TA"):
            return "A"
        return None
