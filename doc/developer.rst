=====================
Developer Information
=====================


Components
----------

This package consists of two components:

* a grammar for HGVS, based on `Parsley`_, wrapped in a simple HGVSParser
  class
* simple classes that model elements of HGVS concepts and provides
  formatting

The current version of the parser requires all components.

Gene names are not accepted, either in lieu of or in addition to the
sequence name (the Recommendations permit both). The recommendation uses
parenthesis to express uncertainty in position, edit, and other details;
this is not supported.


Background
----------

The Human Genome Variation Society, `HGVS`_, produced the `HGVS
Recommendations`_ for the representation of variants in genome,
transcript, and protein sequences.  The following is a *brief* overview of
the standard and especially the subset parsed by this package.


Classes
-------

* hgvs.parser.HGVSParser

* hgvs.edit.DelIns

* hgvs.location.Interval


.. _`Parsley`: https://pypi.python.org/pypi/Parsley
.. _`HGVS`: http://www.hgvs.org/
.. _`HGVS Recommendations`: http://hgvs.org/mutnomen/
