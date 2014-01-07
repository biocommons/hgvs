Overview
========


.. include:: description.rst


Components
----------

This package consists of two components:

* a grammar for HGVS, based on `Parsley`_, wrapped in a simple HGVSParser
  class
* classes that model elements of HGVS concepts

We have made an intentional choice to not implement components of HGVS
that are error prone or ambiguous.  Omitted features include gene names as
a sequence reference (versioned accessions are the only unambiguous
reference) and the use of parentheses to express uncertainty in position,
edit, and other details.


Background
----------

The Human Genome Variation Society, `HGVS`_, produced the `HGVS
Recommendations`_ for the representation of variants in genome,
transcript, and protein sequences.  The following is a *brief* overview of
the standard and especially the subset parsed by this package.



.. _`Parsley`: https://pypi.python.org/pypi/Parsley
.. _`HGVS`: http://www.hgvs.org/
.. _`HGVS Recommendations`: http://hgvs.org/mutnomen/
