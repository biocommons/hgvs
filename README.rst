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

Prelim::

  In [1]: import hgvs.parser
  
  In [2]: hgvs.__version__
  Out[2]: '0.0.4'
  
  In [3]: parser = hgvs.parser.Parser()


Simple stuff::

  In [4]: parser.accn('NM_01234.5')
  Out[4]: 'NM_01234.5'

  In [17]: parser.accn('NM_01234.5.2')
  ParseError: 
  NM_01234.5.2
             ^
  Parse error at line 1, column 11: expected EOF. trail: [digit]

  In [7]: pos = parser.c_interval('12-34_56+78')
  
  In [16]: pos
  Out[16]: Interval(start=12-34, end=56+78)


HGVS Variants::

  In [7]: var = parser.hgvs_variant('NM_01234.5:c.23A>T')
  
  In [8]: var.seqref, var.type, var.posedit
  Out[8]: ('NM_01234.5', 'c', PosEdit(pos=23, edit=A>T))
  
  In [9]: var.posedit.pos
  Out[9]: Interval(start=23, end=23)
  
  In [12]: var.posedit.edit
  Out[12]: DelIns(ref=A, alt=T)
  
  In [13]: var.posedit.edit.ref
  Out[13]: 'A'
  
  In [14]: var.posedit.edit.alt
  Out[14]: 'T'
  
  In [15]: str(var)
  Out[15]: 'NM_01234.5:c.23A>T'
  
  In [18]: var.posedit.edit.alt = 'TT'
  
  In [19]: str(var)
  Out[19]: 'NM_01234.5:c.23delAinsTT'

HGVS Positions::

  In [25]: hp = parser.hgvs_position('NM_01234.5:c.23-12_78+56')
  
  In [26]: hp.pos.start.offset
  Out[26]: -12

  (and surely you get the idea by now)


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
