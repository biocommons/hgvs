.. _installation:
.. _Installing hgvs:

Installing hgvs
!!!!!!!!!!!!!!!


Supported Platforms
@@@@@@@@@@@@@@@@@@@

`hgvs` is developed primarily on Ubuntu systems and has been reported
to work on Mac.  Other platforms and dependency versions are expected
to work but have not been tested. Reports of successful operation on
other platforms (and patches to enable this) are appreciated.
**Python >=3.5 is now required.**



Install Prerequisites
@@@@@@@@@@@@@@@@@@@@@

`hgvs` currently requires PostgreSQL client libraries.  On Ubuntu,
try::

  apt-get install libpq-dev

On a Mac with homebrew::

  brew install postgresql


Use a virtual environment
@@@@@@@@@@@@@@@@@@@@@@@@@

Users are encouraged to use a virtual environment.  The most basic
method for this is::

  $ python3 -m venv venv
  $ source venv/bin/activate

Your shell prompt will change upon activation. 
  
See `this tutorial
<https://realpython.com/python-virtual-environments-a-primer/>`__ for
more information about virtual environments.


Installing hgvs from PyPI (preferred)
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Install hgvs via pip::

  $ pip install hgvs

`hgvs` will install dependencies automatically.



Installing hgvs from source (for developers)
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

For the project at https://github.com/biocommons/hgvs. 

Fetch the source code::

  $ git clone https://github.com/<your github username>/hgvs

Then::

  $ source venv/bin/activate # replace with path to your virtual env
  $ cd hgvs
  $ make develop


.. _seqrepo_install:

Installing SeqRepo (optional)
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

`seqrepo <https://github.com/biocommons/biocommons.seqrepo>`__
provides an easy and efficient mechanism to maintain a local
sequence database.

Install seqrepo::

  $ pip install biocommons.seqrepo

Then, choose a file path that has at least 10GB of space available.
By default, seqrepo will use /usr/local/share/serepo/.  Make that
directory::

  $ mkdir /usr/local/share/seqrepo

Download an instance of the human sequence set::

  $ seqrepo -r /usr/local/share/seqrepo pull

You can skip the -r if you use the default
/usr/local/share/seqrepo/.  This step will take 10-30 minutes, or
more for slow connections.

As with UTA, you tell hgvs to use this feature via an environment
variable::

  $ export HGVS_SEQREPO_DIR=/usr/local/share/seqrepo/20160906


.. _uta_docker_install:
.. _uta_docker:

Local Installation of UTA (optional)
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


The easiest way to install UTA locally is to use the docker image:

  $ docker run -d --name uta_20170117 -p 15032:5432 biocommons/uta:uta_20170117

If you do this, then set::

  $ export UTA_DB_URL=postgresql://anonymous@localhost:15032/uta/uta_20170117

If you don't set this variable, `hgvs` will use the remote uta
database.



Test your installation
@@@@@@@@@@@@@@@@@@@@@@

`hgvs` installs `hgvs-shell`, a command line tool based on
IPython.  It's a convenience utility that imports and initializes
frequently-used components.  Try this::
  
  (default-2.7) snafu$ hgvs-shell
  INFO:root:Starting hgvs-shell 1.0.0a1
  INFO:biocommons.seqrepo:biocommons.seqrepo 0.3.1
  INFO:hgvs.dataproviders.seqfetcher:Using SeqRepo(/usr/local/share/seqrepo/master) sequence fetching
  INFO:hgvs.dataproviders.uta:connected to postgresql://anonymous:anonymous@localhost/uta_dev/uta_20170117...

  In [1]: v = hp.parse_hgvs_variant("NM_033089.6:c.571C>G")

  In [2]: am37.c_to_g(v)
  INFO:biocommons.seqrepo.fastadir.fastadir:Opening for reading: /usr/.../1472015601.985206.fa.bgz
  Out[2]: SequenceVariant(ac=NC_000020.10, type=g, posedit=278801C>G)

  In [3]: am38.c_to_g(v)
  INFO:biocommons.seqrepo.fastadir.fastadir:Opening for reading: /usr/.../1472026864.4364622.fa.bgz
  Out[3]: SequenceVariant(ac=NC_000020.11, type=g, posedit=298157C>G)


Package Versioning
@@@@@@@@@@@@@@@@@@

`hgvs` uses `semantic versioning`_.  For a version `x.y.z`,
incrementing x, y, or z denotes backward-incompatible changes, feature
additions, and bug fixes respectively.

Version numbers for released code come directly from the repository
tag.  Therefore, PyPI version 0.1.2 corresponds exactly to the
repository commit tagged as 0.1.2.

Users (i.e., non-developers) are encouraged to use the PyPI releases
and to specify versions to stay within minor releases for API
stability. For example, a line like::

  hgvs>=1.0,<2

in setup.py or requirements.txt indicates that version 1.0 (any patch
level) is required, and that future 1.x-series releases are
acceptable.

  
