=============================================================================
*hgvs* - Python library to parse, format, validate, and map sequence variants
=============================================================================

This package provides a Python library to facilitate the use of genome,
transcript, and protein variants that are represented using the Human
Genome Variation Society (`mutnomen`_) recommendations.

| **Learn:**  | `Changelog`_
| **Use:** |pypi_badge|  |install_status|
| **Interact:** `Mailing List`_ | `Report an Issue`_
| **Develop:** `Source`_ (status: |build_status|)

Citation:

| **A Python Package for Parsing, Validating, Mapping, and Formatting Sequence Variants Using HGVS Nomenclature.**
| Hart RK, Rico R, Hare E, Garcia J, Westbrook J, Fusaro VA.
| *Bioinformatics*. 2014 Sep 30. `PubMed <http://www.ncbi.nlm.nih.gov/pubmed/25273102>`_ | `Open Access PDF <http://bioinformatics.oxfordjournals.org/content/31/2/268.full.pdf>`_

----

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
  (`VariantMapper <http://pythonhosted.org/hgvs/modules.html#hgvs.variantmapper.VariantMapper>`_ and `Projector
  <http://pythonhosted.org/hgvs/modules.html#hgvs.projector.Projector>`_)
* Reliable handling of regions reference-transcript discrepancy (requires UTA_)
* Tools to validate variants (coming soon)
* Pluggable data providers support alternative sources of transcript mapping
  data
* Extensive automated tests, including those for all variant types and
  "problematic" transcripts

----

Important Notes
---------------

* **You are encouraged to** `browse issues
  <https://bitbucket.org/biocommons/hgvs/issues>`_. Please report any
  issues you find.
* **Use a pip package specification to ensure stay within minor
  releases for API stability.** For example, ``hgvs >=0.3,<0.4``.
* **The default branch is development.** Pulling from default will get
  you a *development* version.  Release versions are determined by
  tags; updating to a specific version (*e.g.,* ``hg up -r 0.3.0``)
  will get you exactly that version as on PyPI.  If you implement a
  new feature, please create an issue first and work in a feature
  branch named like '44-normalization'.


----

An Example
----------

See `Installation instructions
<http://pythonhosted.org/hgvs/using_hgvs.html#installation>`_ if you
have installation troubles.

::

  $ mkvirtualenv hgvs-test
  (hgvs-test)$ pip install --upgrade setuptools
  (hgvs-test)$ pip install hgvs
  (hgvs-test)$ python

  >>> import hgvs.dataproviders.uta
  >>> import hgvs.parser
  >>> import hgvs.variantmapper

  # start with these variants as strings
  >>> hgvs_g, hgvs_c = 'NC_000007.13:g.36561662C>T', 'NM_001637.3:c.1582G>A'

  # parse the genomic variant into a Python structure
  >>> hp = hgvs.parser.Parser()
  >>> var_g = hp.parse_hgvs_variant(hgvs_g)
  >>> var_g
  SequenceVariant(ac=NC_000007.13, type=g, posedit=36561662C>T)

  # SequenceVariants are composed of structured objects, e.g.,
  >>> var_g.posedit.pos.start
  SimplePosition(base=36561662, uncertain=False)

  # format by stringification 
  >>> str(var_g)
  'NC_000007.13:g.36561662C>T'

  # initialize the mapper for GRCh37 with splign-based alignments
  >>> hdp = hgvs.dataproviders.uta.connect()
  >>> evm = hgvs.variantmapper.EasyVariantMapper(hdp,
  ...          primary_assembly='GRCh37', alt_aln_method='splign',
  ...          replace_reference=True)
  
  # identify transcripts that overlap this genomic variant
  >>> transcripts = evm.relevant_transcripts(var_g)
  >>> sorted(transcripts)
  ['NM_001177506.1', 'NM_001177507.1', 'NM_001637.3']

  # map genomic variant to one of these transcripts
  >>> var_c = evm.g_to_c(var_g, 'NM_001637.3')
  >>> var_c
  SequenceVariant(ac=NM_001637.3, type=c, posedit=1582G>A)
  >>> str(var_c)
  'NM_001637.3:c.1582G>A'

  # CDS coordinates use BaseOffsetPosition to support intronic offsets
  >>> var_c.posedit.pos.start
  BaseOffsetPosition(base=1582, offset=0, datum=1, uncertain=False)


There are `more examples in the documentation <http://pythonhosted.org/hgvs/examples.html>`_.

----

Contributing
------------

The hgvs package is intended to be a community project that
facilitates the reliable use of sequence variants.  Code and
documentation contributions are appreciated!  Contributions are more
likely to be quickly incorporated if they:

* are submitted against the default branch head (or close to
  it);
* are made via pull requests;
* are in a named branch, named like <issue#>-meaningful-name;
* are narrowly focused on the bug/feature described by the issue
* have discrete commits with good log messages;
* have tests;
* are formatted code with yapf;
* are generally consistent with the (loose) style of the current code
  with respect to variable naming, etc.


.. _changelog: http://pythonhosted.org/hgvs/changelog.html
.. _documentation: http://pythonhosted.org/hgvs/
.. _invitae: http://invitae.com/
.. _mutnomen: http://www.hgvs.org/mutnomen/
.. _source: https://bitbucket.org/biocommons/hgvs/
.. _uta: http://bitbucket.org/biocommons/uta/
.. _mailing list: https://groups.google.com/forum/#!forum/hgvs-discuss
.. _report an issue: https://bitbucket.org/biocommons/hgvs/issues?status=new&status=open


.. |rtd_badge| image:: https://readthedocs.org/projects/hgvs/badge/?version=latest
  :target: http://hgvs.readthedocs.org/
  :align: middle

.. |pypi_badge| image:: https://badge.fury.io/py/hgvs.png
  :target: https://pypi.python.org/pypi?name=hgvs
  :align: middle

.. |build_status| image:: https://drone.io/bitbucket.org/biocommons/hgvs/status.png
  :target: https://drone.io/bitbucket.org/biocommons/hgvs
  :align: middle 

.. |install_status| image:: https://travis-ci.org/reece/hgvs-integration-test.png?branch=master
  :target: https://travis-ci.org/reece/hgvs-integration-test
  :align: middle

.. http://badge.fury.io/for/py/uta

