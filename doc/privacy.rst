.. _privacy:

Privacy Issues
!!!!!!!!!!!!!!

This page provides details about how the hgvs package works, with a
focus on privacy issues that users may have.  The intent is to provide
users with enough information to assess privacy concerns for their
institutions.


What's not done
@@@@@@@@@@@@@@@

No biologically-relevant data are collected or aggregated from any use
of the hgvs package for any purpose.  Furthermore, variant
manipulation is entirely local. Sequence variants are never sent over
the network.

Someg `hgvs` operations require additional data.  For example, mapping
variants between a genomic reference and a transcript requires
transcript-specific alignment information.  Currently, fetching
addition data requires a network connection.

(We are considering whether and how to provide fully self-contained
installations and do not require network access, but such is not
available at this time.)


.. _dataprovider_queries:

Data Provider Queries
@@@@@@@@@@@@@@@@@@@@@

`hgvs` requires a lot of specialized addition data to validate,
normalize, and map variants.  *All* queries for data are consolidated
into a data provider interface that consists of 11 queries.  The
method signatures, including input arguments, are shown below with a
discussion about privacy consequences.

    fetch_seq(ac, start_i, end_i)
      This method fetches reference sequence in the context of the
      variant and is required in order to validate, normalize, and
      replace variant reference sequences.  By sending accession and
      coordinates, it reveals a specific region of interest (and
      therefore genes and possible clinical conditions).

      The current implementation, which fetches transcripts and
      genomic sequences from UTA, NCBI, and Ensembl, is a measure
      until we complete a comprehensive sequence archive.
    
    data_version(), schema_version()
      Queries for meta data about the data provider.
    
    get_acs_for_protein_seq(seq), get_gene_info(gene), get_tx_exons(tx_ac, alt_ac, alt_aln_method), get_tx_for_gene(gene), get_tx_identity_info(tx_ac), get_tx_info(tx_ac, alt_ac, alt_aln_method), get_tx_mapping_options(tx_ac), get_tx_seq(ac)
      For all of these queries, the inputs are combinations of transcript
      accession, reference accession, gene name. These are likely too
      broad to constitute serious privacy concerns.


Information about current connections
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

The following is an example of the kinds of information available
about a current connection as collected by PostgreSQL.

    ===================  =================================================================================
    ===================  =================================================================================
      datname            uta                                                                              
      usename            anonymous                                                                        
      application_name   hgvs-shell/0.4.0rc2.dev20+n97ead5bf0fed.d20150831                                
      client_addr        162.217.73.242                                                                   
      client_hostname    invitae.static.monkeybrains.net                                                  
      client_port        38318                                                                            
      backend_start      2015-08-31 22:58:26.411654+00                                                    
      query_start        2015-08-31 22:58:30.669956+00                                                    
      state_change       2015-08-31 22:58:30.673533+00                                                    
      waiting            f                                                                                
      state              idle                                                                             
      query              select *                                                                         
                         from tx_exon_aln_v                                                               
                         where tx_ac='NM_170707.3' and alt_ac='NC_000001.10' and alt_aln_method='splign'  
                         order by alt_start_i                                                             
    ===================  =================================================================================

Several of these merit discussion.

    application_name
      Upon connection using the UTA data provider, a string containing the
      name of the python script and hgvs version are passed to the
      postgresql server.  The string typically looks like
      ``hgvs-shell/0.4.0rc2.dev20+n97ead5bf0fed.d20150831``.  Clients may
      override the application_name when calling connect().
    
    client_addr and client_hostname
      The source IP and hostname are available for current
      connections. For most clients, this will mean identifying an
      institution but not specific computers or individuals.
      
    query
      The current or most recently executed query is visible. When
      accessed through the data provider, this field is limited to
      :ref:`dataprovider_queries`.


Historical connection information
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Although we do have historical logs for database connections, they
provide only date, time, and database connection.  Currently, we do
not log queries, although we might choose to periodically log
certain queries for performance monitoring.

