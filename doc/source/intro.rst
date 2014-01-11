Introduction
------------

About HGVS
~~~~~~~~~~

Genome, transcript, and protein sequence variants are typically reported
using the `mutation nomenclature ("mutnomen") recommendations
<http://www.hgvs.org/mutnomen/>`_ provided by the `Human Genome Variation
Society (HGVS) <http://www.hgvs.org/>`_.  If you're reading this, you've
likely already seen HGVS-formatted variants: they usually look as simple
as NM_021960.4:c.740C>T, but can be much more complex.

The mutnomen standard was formulated in an era of specialized sequencing
and cytogenetic analyses, long before the computational analysis of
high-throughput sequencing annotation was envisioned.  Unfortunately, the
complexity of biological phenomena and the breadth of the mutnomen
standard makes it difficult to implement the standard in software, which
in turn makes using the standard in high-throughput analyses difficult.

This package, ``hgvs``, is an easy-to-use Python library for parsing,
representing, formatting, and mapping variants between genome, transcript,
and protein sequences.  The current implementation handles most (but not
all) of the mutnomen standard for precisely defined sequence variants.
The intent is to centralize the subset of HGVS variant manipulation that
is routinely used in modern, high-throughput sequencing analysis.


Components
~~~~~~~~~~

The ``hgvs`` package consists of several interoperable components:

* :class:`hgvs.parser.Parser`, based on Parsley_ and a Parsing Expression Grammar (PEG_).
* :doc:`Classes <modules>` that model HGVS concepts, such as CDS positions with offsets
* formatters that generate HGVS strings from classes
* tools to map variants between genome, transcript, and protein sequences

We have made an intentional choice to not implement some components of
HGVS, especially those that do not refer to variants precisely.  For
instance, our grammar does not permit referring to a variant by gene name
because.



Architecture
~~~~~~~~~~~~

.. sidebar:: Package Architecture in a Nutshell

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




.. _`Parsley`: https://pypi.python.org/pypi/Parsley
.. _`HGVS`: http://www.hgvs.org/
.. _`HGVS Recommendations`: http://hgvs.org/mutnomen/
.. _PEG: http://en.wikipedia.org/wiki/Parsing_expression_grammar
