Development
-----------

This section is intended for developers seeking to extend the hgvs
package.  You should be familiar with the architecture, conventions,
and basic functionality elsewhere in this documentation.


Get Cozy with make
~~~~~~~~~~~~~~~~~~

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

  :upload, upload_docs:
     Upload package to PyPI or docs to pythonhosted.org.




Installation for Development
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~

Yes! We'll be thrilled to have your contributions!

The preferred way to submit a patch is by forking the project on
BitBucket, commiting your changes there, then sending a pull request.

If you have a really worthwhile patch, we'll probably accept a
diff-formatted patch, but that'll make it harder for us and impossible
for you to get credit.
