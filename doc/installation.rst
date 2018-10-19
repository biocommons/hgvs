.. _installation:
.. _Installing hgvs:


Installing hgvs
!!!!!!!!!!!!!!!


Supported Platforms
@@@@@@@@@@@@@@@@@@@

**Operating Systems: Ubuntu, MacOS**

**Python ver. 2.7 required**

Reports of successful operation onother platforms (and patches to enable this) are appreciated.

**Coming soon:** Python 3+s

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


Install Prerequisites
@@@@@@@@@@@@@@@@@@@@@

`hgvs` currently requires PostgreSQL client libraries.  
On Ubuntu::

  apt-get install libpq-dev

On Mac::

  brew install postgresql


Use a virtual environment
@@@@@@@@@@@@@@@@@@@@@@@@@

Users are encouraged to use virtualenv.

On Ubuntu::

  $ sudo apt-get install virtualenvwrapper
  $ mkvirtualenv hgvs-test

Users can use pipenv and conda virtualenv methods as well.


.. note::
   Users (non-developers) refer to "Installing from PyPI".


Installing from PyPI
@@@@@@@@@@@@@@@@@@@@

While in your virtualenv, update `setuptools`::

  $ pip install setuptools --upgrade

Install `hgvs` (will install dependencies automatically)::

  $ pip install hgvs
`Dev Notes: Using conda need to also install Sphinx and sphinx_rtd_theme.`  


Installing from source
@@@@@@@@@@@@@@@@@@@@@@


Fetch the source code::

  $ git clone https://github.com/biocommons/hgvs
  $ cd hgvs
  $ make install


.. _seqrepo_install:

Installing SeqRepo (recommended)
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

`seqrepo <https://github.com/biocommons/biocommons.seqrepo>`__
provides an easy and efficient mechanism to maintain a local
sequence database. Choose a file path that has at least 10GB of space available.

Install seqrepo::

  $ pip install biocommons.seqrepo

Default for `seqrepo` is /usr/local/share/serepo/.  

Execute the following accordingly::

  $ mkdir /usr/local/share/seqrepo
  $ seqrepo -r /usr/local/share/seqrepo pull

This step will take 10-30 minutes, or more for slow connections.

Create a variable for `hgvs`::

  $ export HGVS_SEQREPO_DIR=/usr/local/share/seqrepo/20160906


.. _uta_docker_install:
.. _uta_docker:

Local Installation of UTA (optional)
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


If you do not install `UTA`, `hgvs` will use a remote `UTA` database 

This section is for docker users (preferred `UTA` access)::

  $ docker run -d --name uta_20170117 -p 15032:5432 biocommons/uta:uta_20170117

Set the following variable::

  $ export UTA_DB_URL=postgresql://anonymous@localhost:15032/uta/uta_20170117


Test your installation
@@@@@@@@@@@@@@@@@@@@@@

While in your `hgvs` virtualenv, launch `hgvs-shell`::
  
(hgvs) $ hgvs-shell
Python 2.7.15 | packaged by conda-forge | (default, Oct 12 2018, 14:10:50) 
Type "copyright", "credits" or "license" for more information...
  In [1]: v = hp.parse_hgvs_variant("NM_033089.6:c.571C>G")
  In [2]: v
  Out[2]: SequenceVariant(ac=NM_033089.6, type=c, posedit=571C>G)
  In [3]: am37.c_to_g(v)
  Out[3]: SequenceVariant(ac=NC_000020.10, type=g, posedit=278801C>G)
  In [4]: am38.c_to_g(v)
  Out[4]: SequenceVariant(ac=NC_000020.11, type=g, posedit=298157C>G)



