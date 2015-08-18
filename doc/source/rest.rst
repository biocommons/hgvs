REST interface API reference
~~~~~~~

+--------------------------+----------------------------+------------------------------+-----------------------------------------+
| Endpoint                 | Request arguments          | Respond fields               |Description                              |
+==========================+============================+==============================+=========================================+
| data_version             | (None)                     | data_version                 | UTA data version.                       |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
+--------------------------+----------------------------+------------------------------+-----------------------------------------+
| schema_version           | (None)                     | schema_version               | database schema version                 |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
+--------------------------+----------------------------+------------------------------+-----------------------------------------+
| tx_exons                 | tx_ac                      | tes_exon_set_id              | return transcript exon info for         |
|                          |                            | aes_exon_set_id              | supplied accession.                     |
|                          |                            | tx_ac                        |                                         |
|                          |                            | alt_ac                       |                                         |
|                          |                            | alt_strand                   |                                         |
|                          |                            | alt_aln_method               |                                         |
|                          |                            | ord                          |                                         |
|                          |                            | tx_exon_id                   |                                         |
|                          |                            | alt_exon_id                  |                                         |
|                          |                            | tx_start_i                   |                                         |
|                          |                            | tx_end_i                     |                                         |
|                          |                            | alt_start_i                  |                                         |
|                          |                            | alt_end_i                    |                                         |
|                          |                            | cigar                        |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
+--------------------------+----------------------------+------------------------------+-----------------------------------------+
| tx_info                  | tx_ac                      | alt_aln_method               | return a single transcript info for     |
|                          |                            | cds_start_i                  | supplied accession.                     |
|                          |                            | alt_ac                       |                                         |
|                          |                            | cds_end_i                    |                                         |
|                          |                            | tx_ac                        |                                         |
|                          |                            | hgnc                         |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
+--------------------------+----------------------------+------------------------------+-----------------------------------------+
| sequence                 | ac                         | seq                          | Fetches sequence by accession,          |
|                          |                            |                              | optionally bounded by [start, end) .    |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
+--------------------------+----------------------------+------------------------------+-----------------------------------------+
| tx_for_gene              | gene                       | alt_aln_method               | return transcript info records for      |
|                          |                            | cds_start_i                  | supplied gene, in order of decreasing   |
|                          |                            | alt_ac                       | length.                                 |
|                          |                            | cds_end_i                    |                                         |
|                          |                            | tx_ac                        |                                         |
|                          |                            | hgnc                         |                                         |
|                          |                            |                              |                                         |
+--------------------------+----------------------------+------------------------------+-----------------------------------------+
| tx_for_region            | alt_ac                     | alt_aln_method               | return transcripts that overlap given   |
|                          |                            | start_i                      | region.                                 |
|                          |                            | alt_ac                       |                                         |
|                          |                            | tx_ac                        |                                         |
|                          |                            | end_i                        |                                         |
|                          |                            | alt_strand                   |                                         |
|                          |                            |                              |                                         |
+--------------------------+----------------------------+------------------------------+-----------------------------------------+
| acs_for_protein_seq      | seq                        | ac                           | returns a list of protein accessions    |
|                          |                            |                              | for a given sequence.  The list is      |
|                          |                            |                              | guaranteed to contain at least one      |
|                          |                            |                              | element with the MD5-based accession    |
|                          |                            |                              | (MD5_01234abc...def56789) at the end    |
|                          |                            |                              | of the list.                            |
|                          |                            |                              |                                         |
+--------------------------+----------------------------+------------------------------+-----------------------------------------+
| gene_info                | gene                       | added                        | returns basic information about the     |
|                          |                            | descr                        | gene.                                   |
|                          |                            | summary                      |                                         |
|                          |                            | maploc                       |                                         |
|                          |                            | hgnc                         |                                         |
|                          |                            | aliases                      |                                         |
|                          |                            |                              |                                         |
+--------------------------+----------------------------+------------------------------+-----------------------------------------+
| tx_mapping_options       | tx_ac                      | alt_ac                       | Return all transcript alignment sets    |
|                          |                            | tx_ac                        | for a given transcript accession        |
|                          |                            | alt_aln_method               | (tx_ac); returns empty list if          |
|                          |                            |                              | transcript does not exist.  Use this    |
|                          |                            |                              | method to discovery possible mapping    |
|                          |                            |                              | options supported in the database.      |
|                          |                            |                              |                                         |
+--------------------------+----------------------------+------------------------------+-----------------------------------------+
| tx_identity_info         | tx_ac                      | alt_aln_method               | returns features associated with a      |
|                          |                            | lengths                      | single transcript.                      |
|                          |                            | tx_ac                        |                                         |
|                          |                            | alt_ac                       |                                         |
|                          |                            | cds_end_i                    |                                         |
|                          |                            | cds_start_i                  |                                         |
|                          |                            | hgnc                         |                                         |
+--------------------------+----------------------------+------------------------------+-----------------------------------------+
| similar_transcripts      | tx_ac                      | es_fp_eq                     | Return a list of transcripts that are   |
|                          |                            | cds_es_fp_eq                 | similar to the given transcript, with   |
|                          |                            | cds_eq                       | relevant similarity criteria.           |
|                          |                            | hgnc_eq                      |                                         |
|                          |                            | tx_ac1                       |                                         |
|                          |                            | tx_ac2                       |                                         |
|                          |                            |                              |                                         |
+--------------------------+----------------------------+------------------------------+-----------------------------------------+
| pro_ac_for_tx_ac         | tx_ac                      | pro_ac                       | Return the (single) associated protein  |
|                          |                            |                              | accession for a given transcript        |
|                          |                            |                              | accession, or None if not found.        |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
|                          |                            |                              |                                         |
+--------------------------+----------------------------+------------------------------+-----------------------------------------+



