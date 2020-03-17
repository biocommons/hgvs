Getting Help
!!!!!!!!!!!!

`hgvs` always works and has no bugs. Furthermore, its interface is
so easy to use that a manual is unnecessary. 

Just kidding.

While hgvs is well-tested and has been used by many groups for several
years, bugs, unexpected behaviors, and usages questions occur.
Fortunately, there's now a small community of people who can help.

If you need help, please read the following sources first.  Then, if
you've still got a question, post to one of them.


If you have questions about the `Variation Nomenclature
Recommendations <http://varnomen.hgvs.org/>`_, consider posting your
questions to the `HGVS Facebook page
<https://www.facebook.com/HGVSmutnomen>`_.



hgvs-discuss Mailing List/Group
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

For general questions, the best source of information is the
hgvs-discuss Google Group
(https://groups.google.com/forum/#!forum/hgvs-discuss).  It is
publicly visible, but posting requires joining in order to control
spam.  The mailing list is the preferred way to reach the `hgvs`
package authors.  (Please do NOT send email directly to authors.)


Gitter Channel
@@@@@@@@@@@@@@

We have a new gitter community at https://gitter.im/biocommons/hgvs.
There's not much use yet, but there's a chance that you could get
real-time replies there.


.. _bug-reports:

Bug Reports
@@@@@@@@@@@

If you think you've got a bug, please report it!  Here are a few tips
to make it more likely that you get a useful reply:

* Use the command-line tool `hgvs-shell` that comes with `hgvs` to
  prepare your bug report.  Using `hgvs-shell` makes it easier for you
  to report the bug and make it easier for developers to understand
  it.

* Take the time to prepare a minimal example that demonstrates the
  bug.  You are unlikely to get a reply if you submit a report that
  includes your own wrappers and tooling.

* Include the bug demonstration as text. A screenshot of a bug report
  is not reproducible.

* Include the values of `hgvs.__version__` and `hgvs.hdp.url`, and
  whether you're using seqrepo. (i.e., whether you specified
  HGVS_SEQREPO_DIR)

* hgvs-shell in an upcoming release will provide much of the above
  information for you, as shown below. Please use it.

* Include an explanation of the result you expected and why.

* Report the bug using github, which requires an account.  If you
  don't have an account (and don't want to create one), sending the
  same information to the mailing list is acceptable.

::

  $ hgvs-shell
  
  ############################################################################
  hgvs-shell -- interactive hgvs
  hgvs version: 1.1.3.dev11+ne7b6a1c3ec7a
  data provider url: postgresql://anonymous:anonymous@uta.biocommons.org/uta/uta_20170117
  schema_version: 1.1
  data_version: uta_20170117
  sequences source: remote (bioutils.seqfetcher)
  
  The following variables are defined:
  * hp -- hgvs parser
  * hdp -- hgvs data provider
  * vm -- VariantMapper
  * am37 -- AssemblyMapper, GRCh37
  * am38 -- AssemblyMapper, GRCh38
  * hv -- hgvs Validator
  * hn -- hgvs Normalizer
  
  hgvs_g, hgvs_c, hgvs_p -- sample variants as hgvs strings
  var_g, var_c, var_p -- sample variants, as parsed SequenceVariants
  
  When submitting bug reports, include the version header shown above
  and use these variables/variable names whenever possible.

  In [1]:
  
