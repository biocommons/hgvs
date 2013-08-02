=========================
HGVS Parser and Formatter
=========================

This package, hgvs-utils, implements parsers and formatters for biological
sequence variants using the `HGVS`_ recommendations.  Currently, the
package implements an incomplete -- but very usable -- subset of the
standards.

.. _HGVS: http://www.hgvs.org/mutnomen/


Example
-------

::

  In [1]: import hgvs.parser
  
  In [2]: hp = hgvs.parser.Parser()
  
  In [3]: var= hp.parse('NM_012345.6:c.12+34_98-76delACinsTGTG')
  
  In [4]: var
  Out[4]: Variant(seqref=NM_012345.6, type=c, pos=12+34_98-76, edit=delACinsTGTG)
  
  In [5]: var.pos.start.offset
  Out[5]: 34


Status
------

This library is under development.  The API, grammer, and object model may
change. Please consider pinning any dependencies to a specific version.

* Implemented:

  * simple object model for HGVS variants

  * grammar for genomic and cDNA variants, including intronic variants

  * formatting of variant objects in HGVS format via __str__

  * supports all definite location types, including negative c. and c.*. 

  * supports subst, del, ins, delins, dup, repeat

* Not implemented (partial list):

  * inversions (easy)

  * tranlocations

  * m., n., r., p.

  * compound variants (e.g., NM_01234.5:c.[56A>T];[64C>G])

  * uncertainty (parens and '?')in position, variant, or range

  * gene names
