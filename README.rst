*hgvs* - manipulate biological sequence variants according to Human Genome Variation Society recommendations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

**Important:** biocommons packages (including hgvs) will begin phasing
out support for Python 2.7 on Jan 1, 2019.  `More
<https://groups.google.com/forum/#!topic/hgvs-discuss/iLUzjzoD-28>`__

The *hgvs* package provides a Python library to parse, format,
validate, normalize, and map sequence variants according to `Variation
Nomenclature`_ (aka Human Genome Variation Society) recommendations.

Specifically, the hgvs package focuses on the subset of the HGVS
recommendations that precisely describe sequence-level variation
relevant to the application of high-throughput sequencing to clinical
diagnostics.  The package does not attempt to cover the full scope of
HGVS recommendations. Please refer to `issues
<https://github.com/biocommons/hgvs/issues>`_ for limitations.


+--------------------+--------------------------------------------------------------------+
| **Information**    | | |rtd|   |changelog|  |github_license|  |binder|                  |
|                    | | |gitter|   |group|     |getting_help|                            |
+--------------------+--------------------------------------------------------------------+
| **Latest Release** | |github_tag|   |pypi_rel|                                          |
+--------------------+--------------------------------------------------------------------+
| **Development**    | | |status_rel|  |coveralls|   |hit| (pip install)                  |
| (master branch)    | | |issues|  |github_open_pr|                                       |
|                    | | |github_stars|  |github_forks|   |github_contrib|                |
+--------------------+--------------------------------------------------------------------+



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

.. code-block:: python

  >>> import hgvs.parser

  # start with these variants as strings
  >>> hgvs_g = 'NC_000007.13:g.36561662C>T'
  >>> hgvs_c = 'NM_001637.3:c.1582G>A'

  # parse the genomic variant into a Python structure
  >>> hp = hgvs.parser.Parser()
  >>> var_g = hp.parse_hgvs_variant(hgvs_g)
  >>> var_g
  SequenceVariant(ac=NC_000007.13, type=g, posedit=36561662C>T, gene=None)

  # SequenceVariants are composed of structured objects, e.g.,
  >>> var_g.posedit.pos.start
  SimplePosition(base=36561662, uncertain=False)

  # format by stringification 
  >>> str(var_g)
  'NC_000007.13:g.36561662C>T'


Projecting ("Mapping") variants between aligned genome and transcript sequences
###############################################################################

`hgvs` provides tools to project variants between genome, transcript,
and protein sequences.  Non-coding and intronic variants are
supported.  Alignment data come from the `Universal Transcript Archive
(UTA) <https://github.com/biocommons/uta/>`__.

.. code-block:: python

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
  SequenceVariant(ac=NM_001637.3, type=c, posedit=1582G>A, gene=None)
  >>> str(var_c)
  'NM_001637.3:c.1582G>A'

  # CDS coordinates use BaseOffsetPosition to support intronic offsets
  >>> var_c.posedit.pos.start
  BaseOffsetPosition(base=1582, offset=0, datum=Datum.CDS_START, uncertain=False)


Translating coding variants to protein sequences
################################################

Coding variants may be translated to their protein consequences.  HGVS
uses the same pairing of transcript and protein accessions as seen in
NCBI and Ensembl.

.. code-block:: python

   # translate var_c to its protein consequence
   # The object structure of protein variants is nearly identical to
   # that of nucleic acid variants and is converted to a string form
   # by stringification. Per HGVS recommendations, inferred consequences
   # must have parentheses to indicate uncertainty.
   >>> var_p = am.c_to_p(var_c)
   >>> var_p
   SequenceVariant(ac=NP_001628.1, type=p, posedit=(Gly528Arg), gene=None)
   >>> str(var_p)
   'NP_001628.1:p.(Gly528Arg)'

   # setting uncertain to False removes the parentheses on the
   # stringified form
   >>> var_p.posedit.uncertain = False
   >>> str(var_p)
   'NP_001628.1:p.Gly528Arg'

   # formatting can be customized, e.g., use 1 letter amino acids to
   # format a specific variant
   # (configuration may also be set globally)
   >>> var_p.format(conf={"p_3_letter": False})
   'NP_001628.1:p.G528R'



Normalizing variants
####################

Some variants have multiple representations due to instrinsic
biological ambiguity (e.g., inserting a G in a poly-G run) or due to
misunderstanding HGVS recommendations.  Normalization rewrites certain
veriants into a single representation.

.. code-block:: python

  # rewrite ins as dup (depends on sequence context)
  >>> import hgvs.normalizer
  >>> hn = hgvs.normalizer.Normalizer(hdp)
  >>> hn.normalize(hp.parse_hgvs_variant('NM_001166478.1:c.35_36insT'))
  SequenceVariant(ac=NM_001166478.1, type=c, posedit=35dup, gene=None)

  # during mapping, variants are normalized (by default)
  >>> c1 = hp.parse_hgvs_variant('NM_001166478.1:c.31del')
  >>> c1
  SequenceVariant(ac=NM_001166478.1, type=c, posedit=31del, gene=None)
  >>> c1n = hn.normalize(c1)
  >>> c1n
  SequenceVariant(ac=NM_001166478.1, type=c, posedit=35del, gene=None)
  >>> g = am.c_to_g(c1)
  >>> g
  SequenceVariant(ac=NC_000006.11, type=g, posedit=49917127del, gene=None)
  >>> c2 = am.g_to_c(g, c1.ac)
  >>> c2
  SequenceVariant(ac=NM_001166478.1, type=c, posedit=35del, gene=None)


There are `more examples in the documentation
<http://hgvs.readthedocs.org/en/stable/examples.html>`_.


Citing hgvs (the package)
@@@@@@@@@@@@@@@@@@@@@@@@@

| **hgvs: A Python package for manipulating sequence variants using HGVS nomenclature: 2018 Update.**
| Wang M, Callenberg KM, Dalgleish R, Fedtsov A, Fox N, Freeman PJ, Jacobs KB, Kaleta P, McMurry AJ, Prlić A, Rajaraman V, Hart RK
| Human Mutation. 2018 `Pubmed <https://www.ncbi.nlm.nih.gov/pubmed/30129167>`__ | `Open Access PDF <https://doi.org/10.1002/humu.23615>`__

| **A Python Package for Parsing, Validating, Mapping, and Formatting Sequence Variants Using HGVS Nomenclature.**
| Hart RK, Rico R, Hare E, Garcia J, Westbrook J, Fusaro VA.
| *Bioinformatics*. 2014 Sep 30. `PubMed <http://www.ncbi.nlm.nih.gov/pubmed/25273102>`__ | `Open Access PDF <http://bioinformatics.oxfordjournals.org/content/31/2/268.full.pdf>`__


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

.. |getting_help| image:: https://img.shields.io/badge/!-help%20me-red.svg
   :target: https://hgvs.readthedocs.io/en/stable/getting_help.html

.. |rtd| image:: https://img.shields.io/badge/docs-readthedocs-green.svg
   :target: http://hgvs.readthedocs.io/

.. |changelog| image:: https://img.shields.io/badge/docs-changelog-green.svg
   :target: https://hgvs.readthedocs.io/en/stable/changelog/

.. |github_license| image:: https://img.shields.io/github/license/biocommons/hgvs.svg
   :alt: GitHub license
   :target: https://github.com/biocommons/hgvs/blob/master/LICENSE)

.. |group| image:: https://img.shields.io/badge/group-hgvs%20discuss-green.svg
   :alt: Mailing list
   :target: https://groups.google.com/forum/#!forum/hgvs-discuss

.. |gitter| image:: https://img.shields.io/badge/chat-gitter-green.svg
   :alt: Join the chat at https://gitter.im/biocommons/hgvs
   :target: https://gitter.im/biocommons/hgvs?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge


.. |github_tag| image:: https://img.shields.io/github/tag/biocommons/hgvs.svg
   :alt: GitHub tag
   :target: https://github.com/biocommons/hgvs

.. |pypi_rel| image:: https://img.shields.io/pypi/v/hgvs.svg
   :target: https://pypi.org/project/hgvs/


.. |status_rel| image:: https://img.shields.io/travis/biocommons/hgvs/master.svg
   :target: https://travis-ci.org/biocommons/hgvs?branch=master

.. |coveralls| image:: https://img.shields.io/coveralls/github/biocommons/hgvs.svg
   :target: https://coveralls.io/github/biocommons/hgvs

.. |issues| image:: https://img.shields.io/github/issues-raw/biocommons/hgvs.svg
   :alt: issues
   :target: https://github.com/biocommons/hgvs/issues

.. |github_open_pr| image:: https://img.shields.io/github/issues-pr/biocommons/hgvs.svg
   :alt: GitHub Open Pull Requests
   :target: https://github.com/biocommons/hgvs/pull/

.. |github_stars| image:: https://img.shields.io/github/stars/biocommons/hgvs.svg?style=social&label=Stars
   :alt: GitHub stars
   :target: https://github.com/biocommons/hgvs/stargazers

.. |github_forks| image:: https://img.shields.io/github/forks/biocommons/hgvs.svg?style=social&label=Forks
   :alt: GitHub forks
   :target: https://github.com/biocommons/hgvs/network

.. |github_contrib| image:: https://img.shields.io/github/contributors/biocommons/hgvs.svg
   :alt: GitHub license
   :target: https://github.com/biocommons/hgvs/graphs/contributors/

.. |install_status| image:: https://travis-ci.org/reece/hgvs-integration-test.png?branch=master
   :target: https://travis-ci.org/reece/hgvs-integration-test

.. |binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/biocommons/hgvs/master?filepath=examples

.. |hit| image:: https://travis-ci.org/biocommons/hgvs-installation-test.svg?branch=master
    :target: https://travis-ci.org/biocommons/hgvs-installation-test	    
