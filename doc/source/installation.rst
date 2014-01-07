.. _installation.rst

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
