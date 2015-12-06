Development
!!!!!!!!!!!

This section is intended for developers seeking to extend the hgvs
package.  You should be familiar with the architecture, conventions,
and basic functionality elsewhere in this documentation.



Get Cozy with make
@@@@@@@@@@@@@@@@@@

The hgvs package includes a GNU Makefile that aids nearly all
developer tasks.  It subsumes much of the functionality in setup.py.
While using the Makefile isn't required to develop, it is the official
way to invoke tools, tests, and other development features. Type
`make` for hgvs-specific help.

Some of the key targets are:

  :develop:
     Prepare the directory for local development.

  :install:
     Install hgvs (as with python setup.py install).

  :test:
     Run the default test suite (~4 minutes).

  :test-quick:
     Run the quick test suite (~35s) of most functionality.

  :clean, cleaner, cleanest:
     Remove extraneous files, leaving a directory in various states of
     tidiness.

  :docs:
     Make the sphinx docs in doc/build/html/.

  :upload:
     Upload package to PyPI




Installation for Development
@@@@@@@@@@@@@@@@@@@@@@@@@@@@

::

  $ hg clone ssh://hg@bitbucket.org/biocommons/hgvs
  $ cd hgvs
  $ make develop


Variables
@@@@@@@@@

The following code variable conventions are used for most of the ``hgvs``
code base.  They should be considered aspirations rather than reality or
policy.  Understanding these conventions will help uses and developers
understand the code.

.. note:: A note on variable suffixes
  If a particular variant type is expected, a suffix is often added to
  variable names. |eg| ``var_c`` in a function argument list signifies
  that a SequenceVariant object with type='c' is expected.

:hgvs*: a string representing an HGVS variant name.  

:var*: a :class:`hgvs.variant.SequenceVariant` object

:pos: 

:posedit: 

:hgvs_position:


.. |eg| replace:: *e.g.,*
.. |ie| replace:: *i.e.,*



Submitting Patches
@@@@@@@@@@@@@@@@@@

Yes! We'll be thrilled to have your contributions!

The preferred way to submit a patch is by forking the project on
BitBucket, commiting your changes there, then sending a pull request.

If you have a really worthwhile patch, we'll probably accept a
diff-formatted patch, but that'll make it harder for us and impossible
for you to get credit.


Developing and Contributing to HGVS
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

* Fork the project at https://bitbucket.org/biocommons/hgvs/

* Clone the project locally with:

    $ hg clone https://bitbucket.org/<your_username>/hgvs

* Create a virtualenv

    $ mkvirtualenv hgvs

* Prepare your environment

    $ make develop

(The Makefile in hgvs wraps functionality in setup.py, and also
provides many useful utilitarian rules. Type ``make`` to see a list of
targets.)

* Code away, then commit and push

    $ hg commit -m 'fixes #141: implements Formatter class'

    $ hg push

* If you'd like to contribute back, submit a pull request on the hgvs
  web site.



Using a local/alternative UTA instance
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

* Install UTA from a PostgreSQL as described at in the UTA_ README.

* Specify an alternate UTA instance.

The easiest way to use a UTA instance other than the default is by
setting UTA_DB_URL.  The format is
``postgresql://<user>:<pass>@<host>/<db>/<schema>``. For example:

   ``postgresql://anonymous:anonymous@uta.biocommons.org/uta/uta_20140210``
  
explicitly selects the public database, and 

   ``postgresql://localhost/uta/uta_20140210``
 
selects a local instance.  Developers can test connectivity like this:

   ``$ UTA_DB_URL=postgresql://localhost/uta/uta_20140210 make test-quick``

See hgvs/dataproviders/uta.py for current UTA database URLs.

.. _uta: http://bitbucket.org/biocommons/uta/


Release Process
@@@@@@@@@@@@@@@

``hgvs`` uses a home-grown tool, ``clogger``, to generate change logs.
This section documents the process.  (Clogger will be released at some
point, but it is currently really only executable by Reece.)

``clogger``\'s primary goal is to propose a preliminary changelog
based on commit messages between specified release tags.  That
``.clog`` file is a simple format like this::

    clog format: 1; -*-outline-*-
    * 0.4.1 (2015-09-14)
    Changes since 0.4.0 (2015-09-09).
    ** Bug Fixes
    *** fixes #274, #275: initialize normalizer with same alt_aln_method as EasyVariantMapper [43e174d6f8af]
    *** fixes #276: raise error when user attempts to map to/from c. with non-coding transcript [3f7b659f4f02]

``.clog`` files should be edited for readability during the release
process and committed to the repo (in ``hgvs/doc/changelog/``).

A Makefile in the same directory generates an ``.rst`` file from the
``.clog``. This file must also be committed to the repo.  This file
becomes the release changelog.

Finally, releases are bundled by major.minor versions in a file like
``0.4.rst`` (no patch level). That file must be edited to include the
new release and committed to the repo.


Specific Example -- 0.4.3 release
#################################

The 0.4.x branch has two recent changes for the 0.4.3 release.  Here's
how the release was prepared::

  hg up 0.4.x
  hg tag 0.4.3cl

  cd doc/changelog
  make 0.4.3cl.clog
  mv 0.4.3cl.clog 0.4.3.clog
  #edit 0.4.3.clog for readability
  make 0.4.3.rst
  #edit 0.4.rst to add 0.4.3 to index

``cd ../..`` (hgvs top-level), then ``hg status`` should now look like::

  M doc/changelog/0.4.rst
  A doc/changelog/0.4.3.clog
  A doc/changelog/0.4.3.rst

Check your work. Type ``make docs``, then view ``build/sphinx/html/changelog/0.4.3.html``.

Now we're ready to finish up::

  hg tag --remove 0.4.3cl
  hg com -m 'added docs for 0.4.3 release'
  hg tag 0.4.3
  hg push
  make upload # (builds distribution and uploads to pypi)
