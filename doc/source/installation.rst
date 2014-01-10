.. _installation.rst

Installation
------------

The following instructions were tested on Ubuntu 13.10 (Python 2.7.5+),
Ubuntu 13.04 (Python 2.7.3), and Ubuntu 12.04 (Python 2.7.3).

Optionally, first build a virtualenv::

  $ sudo apt-get install python2.7 python2.7-dev libpq-dev mercurial virtualenvwrapper
  $ mkvirtualenv hgvs-test

Install via pip::

  (hgvs-test)$ pip install setuptools --upgrade
  (hgvs-test)$ pip install hgvs

Alternatively, test and install from source::

  (hgvs-test)$ hg clone ssh://hg@bitbucket.org/invitae/hgvs
  (hgvs-test)$ cd hgvs
  (hgvs-test)$ make develop
  (hgvs-test)$ make test
  (hgvs-test)$ make install

Other platforms and dependency versions are expected to work but are not
tested.

Test your setup like this::

  (hgvs-test)$ python
  Python 2.7.5+ (default, Sep 19 2013, 13:48:49) 
  [GCC 4.8.1] on linux2
  Type "help", "copyright", "credits" or "license" for more information.
  >>> import hgvs.parser
  >>> hgvs.parser.Parser().parse_hgvs_variant("NM_01234.5:c.12+3A>T")
  SequenceVariant(ac=NM_01234.5, type=c, posedit=12+3A>T)
