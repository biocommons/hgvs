.. _quick_start.rst:

Quick Start
!!!!!!!!!!!

This tutorial provides a comprehensive example of how to use the HGVS
package.  Specifically, we'll:

* install hgvs
* parse a genomic variant
* project the genomic variant to all transcripts
* infer the amino acid changes for coding transcripts

We'll use `rs397509113
<https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=397509113>`_
in BRCA1. This variant is coincident with an exon in 3 coding
transcripts, an intron in 2 other coding transcripts, and a non-coding
transcript.

========================  ================================= =========================
transcript (c.)           protein (p.)                      comment
========================  ================================= =========================
NM_007294.3:c.3844del     NP_009225.1:p.(Glu1282AsnfsTer25) 
NM_007297.3:c.3703del     NP_009228.2:p.(Glu1235AsnfsTer25) 
NM_007300.3:c.3844del     NP_009231.2:p.(Glu1282AsnfsTer25) 
NM_007298.3:c.788-655del  NP_009229.2:p.?                   intronic variant
NM_007299.3:c.788-655del  NP_009230.2:p.?                   intronic variant
NR_027676.1:n.3980del     non-coding                        non-coding transcript
========================  ================================= =========================


Install hgvs
@@@@@@@@@@@@

For this demo, you'll obviously need `hgvs`.  In a reasonably modern
environment, the following should suffice::

  $ pip install hgvs

More detailed installation instructions are in :doc:`installation`.


Start hgvs-shell
@@@@@@@@@@@@@@@@

The :mod:`hgvs` package includes an executable called ``hgvs-shell``,
which sets up :mod:`hgvs` for you.  On the command line, type::

  $ hgvs-shell

This is approximately the same thing as::

  $ IPython
  >>> from hgvs.easy import *

:mod:`hgvs.easy` connects to data sources and initializes commonly used
objects that provide most functionality.  

.. note:: Variant validation, normalization, and projection require
	  access to external data, specifically exon structures,
	  transcript alignments, and protein accessions.  Right now,
	  the only source of this data is via the UTA sister projects.
	  When you import :mod:`hgvs.easy`, you will connect to
	  publicly available data sources.  If you want more
	  information on the architecture of :mod:`hgvs` and UTA, see
	  :doc:`intro`.  See :doc:`installation` for information about
	  installing data sources locally for speed and privacy.


Parse the genomic variant
@@@@@@@@@@@@@@@@@@@@@@@@@

In the ``hgvs-shell``, do::

  >>> var_g = parse("NC_000017.11:g.43091687delC")

.. note:: *All* functionality in :mod:`hgvs` is provided by Python
	  classes.  :mod:`hgvs.easy` exposes common methods with
	  functional forms also, which are used in this quick start
	  guide.  For example, ``parse(...)`` above actually calls
	  ```parser.parse(...)``, where ``parser`` is an instance of
	  the :class:`hgvs.parser.Parser` class.

Parsing a variant results in objects that represent the variant. A
SequenceVariant object is comprised of an accession (``ac``), an HGVS
sequence ``type`` (c,g,m,n,r,p), and 0 or more specific sequence
changes (``posedit`` -- a POSition and EDIt).::

   >>> var_g
   SequenceVariant(ac=NC_000017.11, type=g, posedit=43091687del, gene=None)

The ``posedit`` is itself an object of the :class:`hgvs.posedit.PosEdit` class::

  >>> var_g.posedit
  PosEdit(pos=43091687, edit=del, uncertain=False)

The ``pos`` (position) and ``edit`` attributes are also objects that
can represent intervals and more complex edit operations like indels.
The ``uncertain`` flag enables representation of HGVS uncertainty
(typically with parentheses around the uncertain
component). "stringifying" a variant regenerates an HGVS variant::

  >>> str(var_g)
  'NC_000017.11:g.43091687del'

  >>> "This is a variant: {v}".format(v=var_g)
  'This is a variant: NC_000017.11:g.43091687del'

And, in Python 3, stringification works in f-strings, like so::

  > >> f"{var_g}"
  'NC_000017.11:g.43091687del'


Validating and Normalizing Variants
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

`hgvs` provides functionality to validate and normalize variants::

  >>> normalize(var_g)
  SequenceVariant(ac=NC_000017.11, type=g, posedit=43091688del, gene=None)

  >>> validate(var_g)
  True


Projecting variants between sequences
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

When two sequences have alignments available in , a variant may be
"projected" from one sequence to the other.  ``hgvs`` supports
projecting variants

  * from g to c, n
  * from c to g, n, p
  * from n to c, g

The :class:`hgvs.assemblymapper.AssemblyMapper` class provides a
high-level interface to variant projection. :mod:`hgvs.easy`
initializes AssemblyMapper instances for GRCh37 and GRCh37 as ``am37``
and ``am38`` respectively. For example::

  >>> transcripts = am38.relevant_transcripts(var_g)
  >>> sorted(transcripts)
  ['NM_007294.3', 'NM_007297.3', 'NM_007298.3', 'NM_007299.3', 'NM_007300.3', 'NR_027676.1']

We can now project the genomic variant, ``var_g``, to each of these
transcripts using the ``g_to_t`` function, and the transcript variant
to a protein sequnce using the ``t_to_p`` function.

  >>> for ac in get_relevant_transcripts(var_g):
  ...     var_t = g_to_t(var_g, ac)
  ...     var_p = t_to_p(var_t)
  ...     print("-> " + str(var_t) + " (" + str(var_p) + ") ")
  ...
  -> NM_007294.3:c.3844del (NP_009225.1:p.(Glu1282AsnfsTer25))
  -> NM_007297.3:c.3703del (NP_009228.2:p.(Glu1235AsnfsTer25))
  -> NM_007298.3:c.788-655del (NP_009229.2:p.?)
  -> NM_007299.3:c.788-655del (NP_009230.2:p.?)
  -> NM_007300.3:c.3844del (NP_009231.2:p.(Glu1282AsnfsTer25))
  -> NR_027676.1:n.3980del (non-coding)

In ``hgvs``, the ``t`` type can be either ``c`` or ``n``.  Only
variants on coding sequences (``c.``) can be projected to a protein
sequence.  As a special case, ``t_to_p`` returns "non-coding" when the
input variant is on a non-coding sequence.
