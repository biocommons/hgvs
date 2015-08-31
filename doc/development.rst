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
.........

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
