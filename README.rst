====================================================================
HGVS -- Tools to Parse, Format, and Map Biological Sequence Variants
====================================================================

|build_status| | `Source <https://bitbucket.org/invitae/hgvs>`_ | `Documentation <http://pythonhosted.org/hgvs/>`_ | `Discuss <https://groups.google.com/forum/#!forum/hgvs-discuss>`_ | `Issues <https://bitbucket.org/invitae/hgvs/issues?status=new&status=open>`_

This package provides a Python library to facilitate the use of genome,
transcript, and protein variants that are represented using the Human
Genome Variation Society (`HGVS`_) recommendations.

Features
--------

* `A formal grammar <http://pythonhosted.org/hgvs/grammar.html>`_ for HGVS variant names
* `Classes <http://pythonhosted.org/hgvs/modules.html>`_ that model HGVS
  concepts such as `Interval
  <http://pythonhosted.org/hgvs/modules.html#hgvs.location.Interval>`_,
  intronic offsets (in `BaseOffsetPosition
  <http://pythonhosted.org/hgvs/modules.html#hgvs.location.BaseOffsetPosition>`_),
  frameshifts, uncertain positions, and types of variation (`hgvs.edit
  <http://pythonhosted.org/hgvs/modules.html#module-hgvs.edit>`_)
* Formatters that generate HGVS strings from internal representations
* Tools to map variants between genome, transcript, and protein sequences
  (`HGVSMapper <http://pythonhosted.org/hgvs/modules.html#hgvs.hgvsmapper.HGVSMapper>`_ and `Projector
  <http://pythonhosted.org/hgvs/modules.html#hgvs.projector.Projector>`_)
* Reliable handling of regions reference-transcript discrepancy (requires UTA_)
* Tools to validate variants (coming soon)
* Support for alternative sources of reference and transcript mapping
  information (via BDI_)
* Extensive automated tests


An Example
----------
::

  $ ipython
  In [1]: import hgvs.parser
  In [2]: hp = hgvs.parser.Parser()
  In [3]: hgvs_g = 'NC_000007.13:g.36561662C>T'
  In [4]: var_g = hp.parse_hgvs_variant(hgvs_g)
  In [5]: var_g
  Out[5]: Variant(ac=NC_000007.13, type=g, posedit=36561662C>T)

  In [6]: import bdi.sources.uta0_pg, hgvs.hgvsmapper
  In [7]: bdi = bdi.sources.uta0_pg.UTA0()
  In [8]: hm = hgvs.hgvsmapper.HGVSMapper( bdi, cache_transcripts=True )
  In [9]: var_c = hm.hgvsg_to_hgvsc( var_g, 'NM_001637.3' )
  In [10]: var_c
  Out[10]: Variant(ac=NM_001637.3, type=c, posedit=1582G>A)

There are `more examples in the documentation <http://pythonhosted.org/hgvs/examples.html>`_.


.. _HGVS: http://www.hgvs.org/mutnomen/
.. _UTA: http://bitbucket.org/invitae/uta
.. _BDI: http://bitbucket.org/invitae/bdi
.. _Invitae: http://invitae.com/


.. |build_status| image:: https://drone.io/bitbucket.org/invitae/hgvs/status.png
  :target: https://drone.io/bitbucket.org/invitae/hgvs
  :align: middle
