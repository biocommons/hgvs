.. _installation.rst

Installation
------------

The following instructions were tested on Ubuntu 13.10 (Python 2.7.5+),
Ubuntu 13.04 (Python 2.7.3), and Ubuntu 12.04 (Python 2.7.3).  Other
platforms and dependency versions are expected to work but are not tested.

Optionally, first build a virtualenv::

  $ sudo apt-get install python2.7 python2.7-dev libpq-dev mercurial virtualenvwrapper
  $ mkvirtualenv hgvs-test

``mkvirtualenv`` will automatically activate your virtualenv and usually
change the prompt to indicate this.

Ensure you have a current ``setuptools`` package::

  $ pip install setuptools --upgrade

You're now ready to install hgvs via pip::

  $ pip install hgvs

``hgvs`` will install dependencies automatically.



Development
-----------

Alternatively, test and install from source::

  $ hg clone ssh://hg@bitbucket.org/invitae/hgvs
  $ cd hgvs
  $ make develop
  $ make test
  $ make install


Test your installation
----------------------

Test your setup like this::

  (hgvs-test)$ python
  Python 2.7.5+ (default, Sep 19 2013, 13:48:49) 
  [GCC 4.8.1] on linux2
  Type "help", "copyright", "credits" or "license" for more information.
  >>> import hgvs.parser
  >>> hgvs.parser.Parser().parse_hgvs_variant("NM_01234.5:c.12+3A>T")
  SequenceVariant(ac=NM_01234.5, type=c, posedit=12+3A>T)
