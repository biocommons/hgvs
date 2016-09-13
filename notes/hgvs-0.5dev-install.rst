Installing hgvs dev
!!!!!!!!!!!!!!!!!!!

Thanks for installing the new version of ``hgvs``. Your help is
invaluable.

You will be installing a dev release.  It is not intended for
production use.

These instructions will guide you through installing and trying new
features of hgvs.


Major new features
@@@@@@@@@@@@@@@@@@

* Uses uta_20160908 by default, which includes updated transcripts for
  GRCh37 and alignments for GRCh38, downloaded from NCBI on
  2016/09/07.

* Enables local sequence sources in lieu of remote sources. This
  requires downloading those sequences using a new facility called
  `seqrepo <https://github.com/biocommons/biocommons.seqrepo>`__.


Instructions
@@@@@@@@@@@@

#. Create a virtual environment

   The easiest way to isolate changes to your Python environment to
   use `virtualenv <https://virtualenv.pypa.io/en/stable/>`__.  I use
   `virtualenvwrapper
   <https://virtualenvwrapper.readthedocs.io/en/latest/>`__, and the
   instructions below assume that.

   Create and activate a virtualenv::

     $ mkvirtualenv hgvs-test


#. Optional: Install the new `docker <https://www.docker.com/>`__
   instance

   Ideally, you will install this locally with docker::

     $ docker run -d --name uta_20160908 -p 60908:5432 biocommons/uta:uta_20160908

   If you do this, then set::

     $ export UTA_DB_URL=postgresql://anonymous@localhost:60908/uta/uta_20160908

   If you don't set this variable, ``hgvs`` will use the remote uta
   database.

     
#. Optional: Install seqrepo and local sequence archive

   `seqrepo <https://github.com/biocommons/biocommons.seqrepo>`__
   provides an easy and efficient mechanism to maintain a local
   sequence database.

   Install seqrepo::

     $ pip install biocommons.seqrepo==0.3.0.dev0

   Verify you got the right package::

     $ seqrepo --version
     0.3.0.dev0

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


#. Install hgvs

   





Feedback requested
@@@@@@@@@@@@@@@@@@

* What did you install? docker uta? seqrepo? both? How easy was it to install?

* What tests did you undertake?

* What's your sense of performance relative to 0.4.x series?

* Configuration is currently via environment variables.  How would you
  prefer to configure 1) UTA URL, 2) seqrepo path?


