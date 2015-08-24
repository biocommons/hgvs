REST interface API reference
...............

+--------------------------+----------------------------+------------------------------+----------------------------------------------------------+-----------------------------------------+
| Endpoint                 | Request arguments          | Respond fields               | Example URL                                              | Description                             |
+==========================+============================+==============================+==========================================================+=========================================+
| data_version             | (None)                     | - data_version               | `http://api.biocommons.org/hgvs/v0/data_version`         | UTA data version.                       |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
+--------------------------+----------------------------+------------------------------+----------------------------------------------------------+-----------------------------------------+
| schema_version           | (None)                     | - schema_version             | `http://api.biocommons.org/hgvs/v0/schema_version`       | database schema version                 |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
+--------------------------+----------------------------+------------------------------+----------------------------------------------------------+-----------------------------------------+
| tx_exons                 | - ac                       | - tes_exon_set_id            | `http://api.biocommons.org/hgvs/v0/tx_exons?tx_ac=NM\    | return transcript exon info for         |
|                          |                            | - aes_exon_set_id            | _199425.2&alt_ac=NC_000020.10&alt_aln\                   | supplied accession.                     |
|                          |                            | - tx_ac                      | _method=splign`                                          |                                         |
|                          |                            | - alt_ac                     |                                                          |                                         |
|                          |                            | - alt_strand                 |                                                          |                                         |
|                          |                            | - alt_aln_method             |                                                          |                                         |
|                          |                            | - ord                        |                                                          |                                         |
|                          |                            | - tx_exon_id                 |                                                          |                                         |
|                          |                            | - alt_exon_id                |                                                          |                                         |
|                          |                            | - tx_start_i                 |                                                          |                                         |
|                          |                            | - tx_end_i                   |                                                          |                                         |
|                          |                            | - alt_start_i                |                                                          |                                         |
|                          |                            | - alt_end_i                  |                                                          |                                         |
|                          |                            | - cigar                      |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
+--------------------------+----------------------------+------------------------------+----------------------------------------------------------+-----------------------------------------+
| tx_info                  | - tx_ac                    | - alt_aln_method             | `http://api.biocommons.org/hgvs/v0/tx_info?tx_ac=\       | return a single transcript info for     |
|                          |                            | - cds_start_i                | NM_199425.2&alt_ac=NC_000020.10&alt_aln_\                | supplied accession.                     |
|                          |                            | - alt_ac                     | method=splign`                                           |                                         |
|                          |                            | - cds_end_i                  |                                                          |                                         |
|                          |                            | - tx_ac                      |                                                          |                                         |
|                          |                            | - hgnc                       |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
+--------------------------+----------------------------+------------------------------+----------------------------------------------------------+-----------------------------------------+
| sequence                 | - ac                       | - seq                        | `http://api.biocommons.org/hgvs/v0/sequence?ac=\         | Fetches sequence by accession,          |
|                          |                            |                              | NM_199425.2`                                             | optionally bounded by [start, end) .    |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
+--------------------------+----------------------------+------------------------------+----------------------------------------------------------+-----------------------------------------+
| tx_for_gene              | - gene                     | - alt_aln_method             | `http://api.biocommons.org/hgvs/v0/tx_for_gene?gene=\    | return transcript info records for      |
|                          |                            | - cds_start_i                | VSX1`                                                    | supplied gene, in order of decreasing   |
|                          |                            | - alt_ac                     |                                                          | length.                                 |
|                          |                            | - cds_end_i                  |                                                          |                                         |
|                          |                            | - tx_ac                      |                                                          |                                         |
|                          |                            | - hgnc                       |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
+--------------------------+----------------------------+------------------------------+----------------------------------------------------------+-----------------------------------------+
| tx_for_region            | - alt_ac                   | - alt_aln_method             | `http://api.biocommons.org/hgvs/v0/tx_for_region?alt_\   | return transcripts that overlap given   |
|                          |                            | - start_i                    | ac=NC_000020.10&alt_aln_method=splign&start\             | region.                                 |
|                          |                            | - alt_ac                     | =100000&end=200000`                                      |                                         |
|                          |                            | - tx_ac                      |                                                          |                                         |
|                          |                            | - end_i                      |                                                          |                                         |
|                          |                            | - alt_strand                 |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
+--------------------------+----------------------------+------------------------------+----------------------------------------------------------+-----------------------------------------+
| acs_for_protein_seq      | - seq                      | - ac                         | `http://api.biocommons.org/hgvs/v0/acs_for_protein\      | returns a list of protein accessions    |
|                          |                            |                              | _seq?seq=MTTRGFSCLLLLIREIDLSAKRRI`                       | for a given sequence.  The list is      |
|                          |                            |                              |                                                          | guaranteed to contain at least one      |
|                          |                            |                              |                                                          | element with the MD5-based accession    |
|                          |                            |                              |                                                          | (MD5_01234abc...def56789) at the end    |
|                          |                            |                              |                                                          | of the list.                            |
|                          |                            |                              |                                                          |                                         |
+--------------------------+----------------------------+------------------------------+----------------------------------------------------------+-----------------------------------------+
| gene_info                | - gene                     | - added                      | `http://api.biocommons.org/hgvs/v0/gene_info?gene=VSX1`  | returns basic information about the     |
|                          |                            | - descr                      |                                                          | gene.                                   |
|                          |                            | - summary                    |                                                          |                                         |
|                          |                            | - maploc                     |                                                          |                                         |
|                          |                            | - hgnc                       |                                                          |                                         |
|                          |                            | - aliases                    |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
+--------------------------+----------------------------+------------------------------+----------------------------------------------------------+-----------------------------------------+
| tx_mapping_options       | - tx_ac                    | - alt_ac                     | `http://api.biocommons.org/hgvs/v0/tx_mapping_options?\  | Return all transcript alignment sets    |
|                          |                            | - tx_ac                      | tx_ac=NM_199425.2`                                       | for a given transcript accession        |
|                          |                            | - alt_aln_method             |                                                          | (tx_ac); returns empty list if          |
|                          |                            |                              |                                                          | transcript does not exist.  Use this    |
|                          |                            |                              |                                                          | method to discovery possible mapping    |
|                          |                            |                              |                                                          | options supported in the database.      |
|                          |                            |                              |                                                          |                                         |
+--------------------------+----------------------------+------------------------------+----------------------------------------------------------+-----------------------------------------+
| tx_identity_info         | - tx_ac                    | - alt_aln_method             | `http://api.biocommons.org/hgvs/v0/tx_identity_info?\    | returns features associated with a      |
|                          |                            | - lengths                    | tx_ac=NM_199425.2`                                       | single transcript.                      |
|                          |                            | - tx_ac                      |                                                          |                                         |
|                          |                            | - alt_ac                     |                                                          |                                         |
|                          |                            | - cds_end_i                  |                                                          |                                         |
|                          |                            | - cds_start_i                |                                                          |                                         |
|                          |                            | - hgnc                       |                                                          |                                         |
+--------------------------+----------------------------+------------------------------+----------------------------------------------------------+-----------------------------------------+
| similar_transcripts      | - tx_ac                    | - es_fp_eq                   | `http://api.biocommons.org/hgvs/v0/similar_transcripts?\ | Return a list of transcripts that are   |
|                          |                            | - cds_es_fp_eq               | tx_ac=NM_199425.2`                                       | similar to the given transcript, with   |
|                          |                            | - cds_eq                     |                                                          | relevant similarity criteria.           |
|                          |                            | - hgnc_eq                    |                                                          |                                         |
|                          |                            | - tx_ac1                     |                                                          |                                         |
|                          |                            | - tx_ac2                     |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
+--------------------------+----------------------------+------------------------------+----------------------------------------------------------+-----------------------------------------+
| pro_ac_for_tx_ac         | - tx_ac                    | - pro_ac                     | `http://api.biocommons.org/hgvs/v0/pro_ac_for_tx_ac?\    | Return the (single) associated protein  |
|                          |                            |                              | tx_ac=NM_199425.2`                                       | accession for a given transcript        |
|                          |                            |                              |                                                          | accession, or None if not found.        |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
|                          |                            |                              |                                                          |                                         |
+--------------------------+----------------------------+------------------------------+----------------------------------------------------------+-----------------------------------------+



