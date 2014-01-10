.. _architecture.rst:

.. sidebar:: HGVS Architecture in a Nutshell

  :hgvs:
     The ``hgvs`` package is most user-visible component that provides
     parsing, formatting, and mapping.  It depends on the ``bdi`` package
     for access to external data.
  
  :bdi:
     ``bdi`` is an abstract interface to data required for hgvs
     operations.  The intent is to enable access to alternative data
     sources (*e.g.,* NCBI, Ensembl, UCSC, LRG). UTA is currently the only
     implementation of this interface.
  
  :uta:
     ``uta`` is an archive of transcripts and alignments. ``bdi`` supports a
     publicly accessible UTA database (in AWS RDS).  A local sqlite version
     is in testing.

Variant mapping and validation requires access to external data,
specifically exon structures, transcript alignments, and protein
accessions.  In order to isolate the hgvs package from the myriad choices
and tradeoffs, these data are provided through an implementation of the
(abstract) `Bioinformatics Data Interface (BDI)
<http://bitbucket.org/invitae/bdi/>`_.

As of Dec 2013, the only available `BDI
<http://bitbucket.org/invitae/bdi>`_ implementation uses the `Universal
Transcript Archive (UTA) <http://bitbucket.org/invitae/uta>`_, a sister
project that provides access to transcripts and genome-transcript
alignments.  `Invitae <http://invitae.com/>`_ provides a public UTA
database instance that is used by default; see the `UTA
<http://bitbucket.org/invitae/uta>`_ page for instructions on installing
your own PostgreSQL or SQLite version.  In the future, other `BDI
<http://bitbucket.org/invitae/bdi>`_ implmentations may be contributed for
other data sources.
