Contributing
!!!!!!!!!!!!

hgvs is intended to be a community project. Contributions of source
code, test cases, and documentation are welcome!

This section describes development conventions and practices for the
hgvs project.  The intention is to help developers get up-to-speed
quickly.


Highlights
@@@@@@@@@@

* Development occurs in the default branch. (Release branches are
  named for the major-minor release, e.g., 0.4.x.)

* Versioning follows `Semantic Versioning`_.

* The code version is determined solely by the hg tags.  This version
  appears in the package name (e.g.,
  ``hgvs-0.4.4-py2.py3-none-any.whl``) and in the version returned by
  ``hgvs.__version__``.  Updating to a specific version (*e.g.,* ``hg
  up -r 0.4.0``) will get you exactly that version.

* All significant development *must* have an associated issue in `hgvs
  issues`_; *create an issue if necessary*. Other changes *may* have
  have an issue. Please develop using a bookmark or branch named for
  the issue, such as ``44-normalization``.

* Pull requests should be narrowly focused around a bug or feature.
  Discrete commits with good log messages facilitate review.  Consider
  collapsing/squashing commits with ``hgvs rebase --collapse ...`` to
  make the PR concise.  Submit PRs against the default branch head (or
  close to it).

* Abide by current `code style`_.  Use ``make reformat`` to reformat all
  code with `yapf`_ prior to submitting a PR.

* Email the `hgvs-discuss`_ mailing list if you have questions.

* Test your code with ``make test`` before you submit a PR.

* Currently, only Python 2.7 is supported. Support for Python 3.5 is
  slated for the next release
  (`#190 <https://github.com/biocommons/hgvs/issues/190/>`__).


A Quick Contribution Example
@@@@@@@@@@@@@@@@@@@@@@@@@@@@

* Fork the project at https://github.com/biocommons/hgvs/

  You will be able to make changes there and then submit your
  contributions for inclusion into the biocommons repo.

.. spacer


* Clone the project locally with

  ::

     $ hg clone https://github.com/<your_username>/hgvs

.. spacer

* Create a virtualenv (recommended)

  ::

     $ mkvirtualenv hgvs

  There are other ways to make python virtual environment. How you do
  this isn't important, but using a virtual environment is good
  practice.

.. spacer

* Prepare your environment

  ::

     $ make develop

  The Makefile in hgvs wraps functionality in setup.py, and also
  provides many useful rules. See `make`_ for more information.

.. spacer

* Make a branch (for significant changes)

  If you expect to change multiple files, please work in a
  branch. Please name the branch like `141-formatter-class`.

.. spacer

* Code, code, code!

  You probably want to test code with::

    $ make test

  See `Local UTA`_ and `make`_ for tips on accellerating testing.

.. spacer

* Reformat code

  This command will reformat the entire package in-place.::

    $ make reformat
    $ hg commit -m 'fixes #141: implements Formatter class'

  Be sure to commit changes afterward!

.. note: Github recognizes "fixes #nnn" and "closes #nnn" as comments
   that close a feature. The preferred use is "fixes" for bugs and
   "closes" for features.

.. spacer

* Commit and push::

  $ make test
  $ hg push

.. spacer

* Submit a pull request at the `hgvs package`_ web site.



.. _Local UTA:

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


.. _make:

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
     Make the sphinx docs in ``docs/build/html/``.




Code Style
@@@@@@@@@@

The package coding style is based roughly on PEP8_, with the following
changes::

  column_limit = 120
  spaces_before_comment = 4
  split_before_named_assigns = True

These code conventions are uniformly enforce by yapf_.  The entire code
base is periodically automatically reformatted for consistency.


Variables
#########

The following code variable conventions are used for most of the `hgvs`
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



Release Process
@@@@@@@@@@@@@@@

`hgvs` uses a home-grown tool, `clogger`, to generate change logs.
This section documents the process.  (Clogger will be released at some
point, but it is currently really only executable by Reece.)

`clogger`\'s primary goal is to propose a preliminary changelog
based on commit messages between specified release tags.  That
``.clog`` file is a simple format like this::

    clog format: 1; -*-outline-*-
    * 0.4.1 (2015-09-14)
    Changes since 0.4.0 (2015-09-09).
    ** Bug Fixes
    *** fixes #274, #275: initialize normalizer with same alt_aln_method as AssemblyMapper [43e174d6f8af]
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

