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

* Uses uta_20160908 and GRCh38 by default. UTA 20160908 includes
  updated transcripts for GRCh37 and alignments for GRCh38, downloaded
  from NCBI on 2016/09/07.  It should be a complete superset of data
  from the previous version, uta_20150903.

* Enables using local sequence sources in lieu of remote sources. This
  requires downloading those sequences using a new facility called
  `seqrepo <https://github.com/biocommons/biocommons.seqrepo>`__.  The
  intended direction is to remove sequences from UTA entirely (it
  contains only transcript sequences currently), and to rely on
  seqrepo (optionally) or remote sequence fetching from NCBI or
  Ensembl.


Instructions
@@@@@@@@@@@@

#. Create a virtual environment

   The easiest way to isolate changes to your Python environment to
   use a virtual environment.

   Preferred: If you use `virtualenvwrapper
   <https://virtualenvwrapper.readthedocs.io/en/latest/>`__::

     $ mkvirtualenv hgvs-test

   Alternative: If you use `virtualenv
   <https://virtualenv.pypa.io/en/stable/>`__::

     $ virtualenv hgvs-test
     $ source hgvs-test/bin/activate

   **With either method, make sure that you're using python 2.**
   (Python 3 support is in the works, and many dependencies have
   already been migrated... it's coming, really!)


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

   The hgvs 0.5.0.devN releases are available in pypi.  Install the lastest like this::

     $ pip install 'hgvs>=0.5.0.dev'
   
#. Try the shell

   hgvs installs ``hgvs-shell``, a command line tool based on IPython.
   It's a convenience utility that imports and initializes
   frequently-used components.  Try this::
     
     (default-2.7) snafu$ hgvs-shell
     INFO:root:Starting hgvs-shell 0.5.0.dev0
     INFO:biocommons.seqrepo:biocommons.seqrepo 0.3.0.dev2
     INFO:hgvs.dataproviders.seqfetcher:Using SeqRepo(/usr/local/share/seqrepo/20160906) sequence fetching
     INFO:hgvs.dataproviders.uta:connected to postgresql://anonymous@localhost:60908/uta/uta_20160908...

     In [1]: v = hp.parse_hgvs_variant("NM_033089.6:c.571C>G")

     In [2]: am37.c_to_g(v)
     INFO:biocommons.seqrepo.fastadir.fastadir:Opening for reading: /usr/.../1472015601.985206.fa.bgz
     Out[2]: SequenceVariant(ac=NC_000020.10, type=g, posedit=278801C>G)

     In [3]: am38.c_to_g(v)
     INFO:biocommons.seqrepo.fastadir.fastadir:Opening for reading: /usr/.../1472026864.4364622.fa.bgz
     Out[3]: SequenceVariant(ac=NC_000020.11, type=g, posedit=298157C>G)


#. Try it on your code

   hgvs 0.5.0 uses GRCh38 **by default**.  You can change that easily
   when invoking AssemblyMapper using the ``assembly_name``
   argument. (NOTE: This command-line argument changed from
   ``primary_assembly`` to ``assembly_name``.  This is the only API
   change, I believe.)


#. IMPORTANT! Send Feedback!

   Please send success and failures to hgvs-discuss@googlegroups.com. 

   In particular, I'd like feedback on the following (at least):

   * For successes and failures, please include OS, Python, and docker
     versions.

   * If you installed the docker UTA instance, was the installation smooth?

   * If you installed seqrepo, was the installation smooth?

   * What tests did you try?  If you had to make significant code
     changes, please describe.

   * Any functional or performance concerns?

   * Configuration is currently via environment variables.  Is that
     acceptable?  If not, how would you prefer to configure 1) UTA
     URL, 2) seqrepo path?

   * All opinions are welcome.  I want to hear any feedback that makes
     ``hgvs`` more useful.
