*hgvs* - manipulate biological sequence variants according to Human Genome Variation Society recommendations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The *hgvs* package provides a Python library to parse, format,
validate, normalize, and map sequence variants according to `Variation
Nomenclature`_ (aka Human Genome Variation Society) recommendations.

+-----------------+-------------------------------------------------------------------------------+
| **Project**     | |rtd_rel| / |gitter| / `Mailing List`_ / `Open Issues`_ / `Apache License`_   |
+-----------------+-------------------------------------------------------------------------------+
| **Release**     | |pypi_rel| / `Change Log`_                                                    |
+-----------------+-------------------------------------------------------------------------------+
| **Development** | `Code`_  / |status_rel|                                                       |
+-----------------+-------------------------------------------------------------------------------+


Features
@@@@@@@@

* Parsing is based on formal grammar.
* An easy-to-use object model that represents
  most variant types (SNVs, indels, dups, inverstions, etc) and
  concepts (intronic offsets, uncertain positions, intervals)
* A variant normalizer that rewrites variants in canoncial forms and
  substitutes reference sequences (if reference and transcript
  sequences differ)
* Formatters that generate HGVS strings from internal representations
* Tools to map variants between genome, transcript, and protein sequences
* Reliable handling of regions genome-transcript discrepancies
* Pluggable data providers support alternative sources of transcript mapping
  data
* Extensive automated tests, including those for all variant types and
  "problematic" transcripts
* Easily installed using remote data sources.  Installation with local
  data sources is straightforward and completely obviates network
  access


Important Notes
@@@@@@@@@@@@@@@

* **You are encouraged to** `browse issues
  <https://github.com/biocommons/hgvs/issues>`_.  All known issues are
  listed there.  Please report any issues you find.
* **Use a pip package specification to stay within minor releases.**
  For example, ``hgvs>=1.1,<1.2``. `hgvs` uses `Semantic Versioning
  <http://semver.org/>`__.


Examples
@@@@@@@@

Installation
#############

By default, `hgvs` uses remote data sources, which makes installation
easy.  

::

  $ mkvirtualenv hgvs-test
  (hgvs-test)$ pip install --upgrade setuptools
  (hgvs-test)$ pip install hgvs
  (hgvs-test)$ python

See `Installation instructions
<http://hgvs.readthedocs.org/en/stable/installation.html>`__ for
details, including instructions for installing `Universal Transcript
Archive (UTA) <https://github.com/biocommons/uta/>`__ and `SeqRepo
<https://github.com/biocommons/biocommons.seqrepo/>`__ locally.


Parsing and Formating
#####################

`hgvs` parses HGVS variants (as strings) into an object model, and can format
object models back into HGVS strings.

::

  >>> import hgvs.parser

  # start with these variants as strings
  >>> hgvs_g = 'NC_000007.13:g.36561662C>T'
  >>> hgvs_c = 'NM_001637.3:c.1582G>A'

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


Projecting ("Mapping") variants between aligned seqeunces
#########################################################

`hgvs` provides tools to project variants between genome, transcript,
and protein sequences.  Non-coding and intronic variants are
supported.  Alignment data come from the `Universal Transcript Archive
(UTA) <https://github.com/biocommons/uta/>`__.

::

  >>> import hgvs.dataproviders.uta
  >>> import hgvs.assemblymapper

  # initialize the mapper for GRCh37 with splign-based alignments
  >>> hdp = hgvs.dataproviders.uta.connect()
  >>> am = hgvs.assemblymapper.AssemblyMapper(hdp,
  ...          assembly_name='GRCh37', alt_aln_method='splign',
  ...          replace_reference=True)
  
  # identify transcripts that overlap this genomic variant
  >>> transcripts = am.relevant_transcripts(var_g)
  >>> sorted(transcripts)
  ['NM_001177506.1', 'NM_001177507.1', 'NM_001637.3']

  # map genomic variant to one of these transcripts
  >>> var_c = am.g_to_c(var_g, 'NM_001637.3')
  >>> var_c
  SequenceVariant(ac=NM_001637.3, type=c, posedit=1582G>A)
  >>> str(var_c)
  'NM_001637.3:c.1582G>A'

  # CDS coordinates use BaseOffsetPosition to support intronic offsets
  >>> var_c.posedit.pos.start
  BaseOffsetPosition(base=1582, offset=0, datum=Datum.CDS_START, uncertain=False)


Normalizing variants
####################

Some variants have multiple representations due to instrinsic
biological ambiguity (e.g., inserting a G in a poly-G run) or due to
misunderstanding HGVS recommendations.  Normalization rewrites certain
veriants into a single representation.

::

  # rewrite ins as dup (depends on sequence context)
  >>> import hgvs.normalizer
  >>> hn = hgvs.normalizer.Normalizer(hdp)
  >>> hn.normalize(hp.parse_hgvs_variant('NM_001166478.1:c.35_36insT'))
  SequenceVariant(ac=NM_001166478.1, type=c, posedit=35dup)

  # during mapping, variants are normalized (by default)
  >>> c1 = hp.parse_hgvs_variant('NM_001166478.1:c.31del')
  >>> c1
  SequenceVariant(ac=NM_001166478.1, type=c, posedit=31del)
  >>> c1n = hn.normalize(c1)
  >>> c1n
  SequenceVariant(ac=NM_001166478.1, type=c, posedit=35del)
  >>> g = am.c_to_g(c1)
  >>> g
  SequenceVariant(ac=NC_000006.11, type=g, posedit=49917127del)
  >>> c2 = am.g_to_c(g, c1.ac)
  >>> c2
  SequenceVariant(ac=NM_001166478.1, type=c, posedit=35del)


There are `more examples in the documentation
<http://hgvs.readthedocs.org/en/stable/examples.html>`_.


Citing hgvs (the package)
@@@@@@@@@@@@@@@@@@@@@@@@@

| **A Python Package for Parsing, Validating, Mapping, and Formatting Sequence Variants Using HGVS Nomenclature.**
| Hart RK, Rico R, Hare E, Garcia J, Westbrook J, Fusaro VA.
| *Bioinformatics*. 2014 Sep 30. `PubMed <http://www.ncbi.nlm.nih.gov/pubmed/25273102>`_ | `Open Access PDF <http://bioinformatics.oxfordjournals.org/content/31/2/268.full.pdf>`_


Contributing
@@@@@@@@@@@@

The hgvs package is intended to be a community project.  Please see
`Contributing
<http://hgvs.readthedocs.org/en/stable/contributing.html>`__ to get
started in submitting source code, tests, or documentation.  Thanks
for getting involved!


See Also
@@@@@@@@

Other packages that manipulate HGVS variants:

* `pyhgvs <https://github.com/counsyl/hgvs>`__
* `Mutalyzer <https://mutalyzer.nl/>`__



.. _docs: http://hgvs.readthedocs.org/
.. _Variation Nomenclature: http://varnomen.hgvs.org/
.. _mailing list: https://groups.google.com/forum/#!forum/hgvs-discuss
.. _Open Issues: https://github.com/biocommons/hgvs/issues
.. _code: https://github.com/biocommons/hgvs
.. _change log: https://hgvs.readthedocs.io/en/stable/changelog/
.. _pypi: https://pypi.python.org/pypi/hgvs
.. _Apache License: https://github.com/biocommons/hgvs/blob/master/LICENSE.txt


.. |rtd_rel| image:: https://readthedocs.org/projects/hgvs/badge/
  :target: http://hgvs.readthedocs.io/
  :align: middle

.. |pypi_rel| image:: https://badge.fury.io/py/hgvs.png
  :target: https://pypi.python.org/pypi?name=hgvs
  :align: middle

.. |status_rel| image:: https://travis-ci.org/biocommons/hgvs.png?branch=master
  :target: https://travis-ci.org/biocommons/hgvs?branch=master
  :align: middle 


.. |install_status| image:: https://travis-ci.org/reece/hgvs-integration-test.png?branch=master
  :target: https://travis-ci.org/reece/hgvs-integration-test
  :align: middle



.. |gitter| image:: https://badges.gitter.im/biocommons/hgvs.svg
   :alt: Join the chat at https://gitter.im/biocommons/hgvs
   :target: https://gitter.im/biocommons/hgvs?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge
