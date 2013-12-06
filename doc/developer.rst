=====================
Developer Information
=====================


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


Classes
-------

* hgvs.parser.HGVSParser

* hgvs.edit
** hgvs.edit.RefAlt

* hgvs.location.Interval


Variable Conventions
--------------------

* var, var_g, var_r, var_c, var_p -- hgvs.variant.Variant instances of unknown or specified/expected types

* hgvs, hgvs_g, hgvs_r, hgvs_c, hgvs_p -- HGVS *strings* of unknown or specified/expected type




.. _`Parsley`: https://pypi.python.org/pypi/Parsley
.. _`HGVS`: http://www.hgvs.org/
.. _`HGVS Recommendations`: http://hgvs.org/mutnomen/
