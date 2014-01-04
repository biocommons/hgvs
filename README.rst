====================================================================
HGVS -- Tools to Parse, Format, and Map Biological Sequence Variants
====================================================================

.. image:: https://drone.io/bitbucket.org/invitae/hgvs/status.png
  :target: https://drone.io/bitbucket.org/invitae/hgvs

This package provides a Python library to facilitate the use of genome,
transcript, and protein variants that are represented using the Human
Genome Variation Society (`HGVS`_) recommendations. ::

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


Features
--------

* Grammar-based parsing of HGVS variants and variant components
* Object model for variants and variant components
* Simple canonicalization of variants
* Mapping of variants between genome, transcript, and protein coordinates (requires `UTA`_; see `Requirements`)
* Extensive automated tests for mapping, including in regions of genome-transcript discrepancies


Status
------

The HGVS package is under active development (Dec 2013).  Although the
grammar, mapping, and testing are fairly mature, API-visible changes are
expected.  We will endeavor to bump the minor version number (in a
major.minor.patch versioning scheme) when we expect code changes are
required.

This packages intends to implement features of `HGVS`_ Nomeclature that are in
widespread use.  Certain common features are not yet implemented, including:

* inversions
* tranlocations
* compound variants

Feedback and bug reports are welcome.


Installation
------------

The following instructions were tested on Ubuntu 13.10 (Python 2.7.5+),
Ubuntu 13.04 (Python 2.7.3), and Ubuntu 12.04 (Python 2.7.3).

First, build a virtualenv::

  $ sudo apt-get install python2.7 python2.7-dev libpq-dev mercurial virtualenvwrapper
  $ mkvirtualenv hgvs-test

Install via pip::

  (hgvs-test)$ pip install setuptools --upgrade
  (hgvs-test)$ pip install hg+ssh://hg@bitbucket.org/invitae/uta
  (hgvs-test)$ pip install hg+ssh://hg@bitbucket.org/invitae/bdi
  (hgvs-test)$ pip install hg+ssh://hg@bitbucket.org/invitae/hgvs

Alternatively, test and install from source::

  (hgvs-test)$ hg clone ssh://hg@bitbucket.org/invitae/hgvs
  (hgvs-test)$ cd hgvs
  (hgvs-test)$ make develop
  (hgvs-test)$ make test
  (hgvs-test)$ make install

Other platforms and dependency versions are expected to work but are not
tested.


Data Dependencies and Architecture
..................................

Variant mapping and validation requires access to external data,
specifically exon structures, transcript alignments, and protein
accessions.  In order to isolate the hgvs package from the myriad choices
and tradeoffs, these data are provided through an implementation of the
(abstract) Bioinformatics Data Interface (`BDI`_).  

As of Dec 2013, the only available BDI implementation uses the Universal
Transcript Archive (`UTA`_), a sister project that provides access to
transcripts and genome-transcript alignments.  `Invitae`_ provides a
public UTA database instance that is used by default; see the `UTA`_
page for instructions on installing your own PostgreSQL or SQLite
version.  In the future, other BDI implmentations may be contributed for
other data sources.

The architecture in a nutshell:

:hgvs:
   The ``hgvs`` package is most user-visible component that provides
   parsing, formatting, and mapping.  It depends on the ``bdi`` package.

:bdi:
   ``bdi`` provider interface for hgvs mapping. Although the intent is to
   enable non-UTA sources eventually, only UTA is currently
   supported. This layer also provides functionality that doesn't belong
   in the ``uta`` package, like fasta file indexing.

:uta:
   ``uta`` is an archive of transcripts and alignments. ``bdi`` supports a
   publicly accessible UTA database (in AWS RDS).  A local sqlite version
   is in testing.




.. _HGVS: http://www.hgvs.org/mutnomen/
.. _UTA: http://bitbucket.org/invitae/uta
.. _BDI: http://bitbucket.org/invitae/bdi
.. _Invitae: http://invitae.com/
