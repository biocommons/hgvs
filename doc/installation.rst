.. _installation:

Installation
!!!!!!!!!!!!


Package Versioning
@@@@@@@@@@@@@@@@@@

The hgvs package uses versions of the form `major.minor.patch`.
Changes in patch level are for bug fixes only and will not have API
changes.  As hgvs is still a nascent project, API changes are possible
even for minor (y) version changes.  (As hgvs matures, will eventually
adopt `semantic versioning`_, in which breaking
API changes will be restricted to major releases. We're not ready for
that yet.)

Version numbers for released code come directly from the repository
tag.  Therefore, PyPI version 0.1.2 corresponds exactly to the
repository commit tagged as 0.1.2.

Users (i.e., non-developers) are encouraged to use the PyPI releases
and to specify versions to stay within minor releases for API
stability. For example, a line like::

  hgvs>=0.4,<0.5

in setup.py or requirements.txt will use the latest patch level
release (with bug fixes!) within the 0.4 series.


Supported Platforms
@@@@@@@@@@@@@@@@@@@

`hgvs` is developed primarily on Ubuntu systems (12.04 through 15.04)
and has been reported to work on Mac.  Other platforms and dependency
versions are expected to work. Reports of successful operation on
other platforms (and patches to enable this) are appreciated.


Install Prerequisites
@@@@@@@@@@@@@@@@@@@@@

``hgvs`` currently requires PostgreSQL client libraries.  (We are
planning to switch to a REST interface and eliminate this dependency
in the 0.5.0 release.)  On Ubuntu, try::

  apt-get install libpq-dev

On a Mac with homebrew::

  brew install postgresql

Users are encouraged to use virtualenv (but this is optional)::

  $ sudo apt-get install python2.7 python2.7-dev libpq-dev mercurial virtualenvwrapper
  $ mkvirtualenv hgvs-test

``mkvirtualenv`` will automatically activate your virtualenv and usually
change the prompt to indicate this.



Installing from PyPI
@@@@@@@@@@@@@@@@@@@@

Ensure you have a current ``setuptools`` package::

  $ pip install setuptools --upgrade

You're now ready to install hgvs via pip::

  $ pip install hgvs

``hgvs`` will install dependencies automatically.



Installing from source
@@@@@@@@@@@@@@@@@@@@@@

.. note::
   Users (non-developers) should prefer the PyPI installation.  There
   is no advantage to installing from source.

Fetch the source code::

  $ git clone https://github.com/biocommons/hgvs

Then::

  $ cd hgvs
  $ make install



Test your installation
@@@@@@@@@@@@@@@@@@@@@@

Test your setup like this::

  (hgvs-test)$ python
  Python 2.7.5+ (default, Sep 19 2013, 13:48:49) 
  [GCC 4.8.1] on linux2
  Type "help", "copyright", "credits" or "license" for more information.
  >>> import hgvs.parser
  >>> hgvs.parser.Parser().parse_hgvs_variant("NM_01234.5:c.12+3A>T")
  SequenceVariant(ac=NM_01234.5, type=c, posedit=12+3A>T)


.. _uta_docker:

Local Installation of UTA
@@@@@@@@@@@@@@@@@@@@@@@@@

The public UTA instance is available without restrictions and is
suitable for many uses.  However, some users may wish to install UTA
locally for performance, isolation, or :ref:`privacy`; see
`instructions for local UTA installation
<https://github.com/biocommons/uta/>`_ for details. 

To override the UTA instance used by the hgvs package, sers must set
the UTA_DB_URL environment variable.  A sample interaction::

  $ docker run --name uta_20160908 -p 60908:5432 biocommons/uta:uta_20160908
  $ export UTA_DB_URL=postgresql://anonymous@localhost:15032/uta/uta_20160908
  $ python -c 'import hgvs.dataproviders.uta; print(hgvs.dataproviders.uta.connect().data_version());'
  uta_20160908

NOTE: The details of UTA_DB_URL string will depend on details of your local
installation.
