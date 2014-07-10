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
because gene names do not precisely define a reference sequence for the
variant.



Architecture
~~~~~~~~~~~~

.. sidebar:: Package Architecture in a Nutshell

  :hgvs:
     The ``hgvs`` package is most user-visible component that provides
     parsing, formatting, and mapping. It requires exon structures and
     sequences to operate, and obtains this data through a class that
     implements the :class:`hgvs.dataproviders.interface`. A concrete
     implementation of this interface uses UTA (below), but users may
     implement their own.
  
  :uta:
     ``uta`` is an archive of transcripts, transcript sequences, and
     transcript-reference sequence alignments.  Invitae provides a
     public UTA instance at ``uta.invitae.com:5432`` (PostgreSQL).

Variant mapping and validation requires access to external data,
specifically exon structures, transcript alignments, and protein
accessions.  In order to isolate the hgvs package from the myriad
choices and tradeoffs, these data are provided through an
implementation of the (abstract) Data Provider Interface
(:class:`hgvs.dataproviders.interface`).


.. _`Parsley`: https://pypi.python.org/pypi/Parsley
.. _`HGVS`: http://www.hgvs.org/
.. _`HGVS Recommendations`: http://hgvs.org/mutnomen/
.. _PEG: http://en.wikipedia.org/wiki/Parsing_expression_grammar
