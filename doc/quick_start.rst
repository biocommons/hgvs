.. _quick_start.rst:

Quick Start
!!!!!!!!!!!

This tutorial provides a comprehensive example of how to use the HGVS
package.  Specifically, we'll:

* install hgvs
* parse a transcript (c.) variant in MCL1 obtained from dbSNP
* project that variant to genomic coordinates (as a g. variant)
* project it back on to another transcript in the same gene
* deduce the amino acid change for that variant

We'll use `rs201430561
<http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=201430561>`_ in
MCL1. This gene has several transcripts, and therefore the genomic variant
NC_000001.10:g.150550916G>A has several distinct transcript
representations:

========================  ==========================
transcript (c.)           protein (p.)
========================  ==========================
NM_001197320.1:c.281C>T   NP_001184249.1:p.Ser94Phe
NM_021960.4:c.740C>T      NP_068779.1:p.Ser247Phe
NM_182763.2:c.688+403C>T  *intronic*
========================  ==========================

This variant was chosen because it has data in dbSNP for comparison and
because it has an intronic variant to spice up the example.


Install ``hgvs``
@@@@@@@@@@@@@@@@

For this demo, you'll need hgvs (of course).  We recommend that you
install IPython as well.  In a reasonably modern environment, the
following should suffice::

  $ pip install --upgrade setuptools
  $ pip install hgvs ipython

More detailed installation instructions are in :doc:`installation`.



Parse the variant
@@@@@@@@@@@@@@@@@

To parse variants, we need to create a an instance of the
:py:class:`hgvs.parser.Parser`.  Since building the grammar is
computationally expensive, you should be only one instance and use it for
all parsing operations.  Start ``ipython``, then do this:

>>> import hgvs.parser
>>> hgvsparser = hgvs.parser.Parser()
>>> var_c1 = hgvsparser.parse_hgvs_variant('NM_001197320.1:c.281C>T')

Parsing a variant results in objects that represent the variant (rather
than, say, Python dictionaries). A SequenceVariant object is comprised of
an accession (``ac``), an HGVS sequence ``type`` (c,g,m,n,r,p), and 0 or
more specific sequence changes (``posedit`` -- a POSition and EDIt).

>>> var_c1
SequenceVariant(ac=NM_001197320.1, type=c, posedit=281C>T)

The ``posedit`` is itself an object:

>>> var_c1.posedit
PosEdit(pos=281, edit=C>T, uncertain=False)

The ``pos`` (position) and ``edit`` attributes are also objects that can
represent intervals and more complex edit operations like indels.  The
``uncertain`` flag enables representation of HGVS uncertainty (typically
with parentheses around the uncertain component).

Finally, "stringifying" a variant regenerates an HGVS variant:

>>> str(var_c1)
'NM_001197320.1:c.281C>T'



Create an VariantMapper instance
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Mapping variants between genomic (g.), transcript (c.), and protein (p.)
sequences is performed by an instance of :py:class:`hgvs.variantmapper.VariantMapper`. As with
the parser, you need only one instance per session.

Variant mapping and validation requires access to external data,
specifically exon structures, transcript alignments, and protein
accessions.  Right now, the only source of this data is via the UTA
UTA sister projects.  (If you want more information on the
architecture of HGVS, UTA, see :doc:`intro`.  However, you don't
really need to understand the architecture to use HGVS.)

First, connect to UTA via :class:``hgvs.dataproviders.uta``:

>>> import hgvs.dataproviders.uta
hdp = hgvs.dataproviders.uta.connect()

By default, you'll connect to the public UTA database instance hosted by
`Invitae <http://invitae.com/>`_.

Then, with that connection, instantiate an VariantMapper:

>>> import hgvs.variantmapper
variantmapper = hgvs.variantmapper.VariantMapper(hdp)

We can use this mapper to transform our transcript variant to a protein variant:

>>> variantmapper.c_to_p(var_c1)
SequenceVariant(ac=NP_001184249.1, type=p, posedit=Ser94Phe)


Map our variant to the genome
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Mapping between sequences is straightforward:

>>> var_g = variantmapper.c_to_g(var_c1,'GRCh37.p10')
>>> var_g
SequenceVariant(ac=NC_000001.10, type=g, posedit=150550916G>A)
>>> str(var_g)
'NC_000001.10:g.150550916G>A'

Notice that this agrees with dbSNP! Also notice that our C>T variant is on
a minus-strand transcript, so the nucleotides are reverse complemented.

Since you'll probably want to access the position, now is a good time to
explore the posedit structure:

First, a posedit consists of a position and an edit.  Positions are
**always** intervals, even if their string representation looks like a
simple integer.  Interval bounds are referred to with ``start`` and
``end`` attributes.  As with edits, they may also be ``uncertain``.
 
>>> var_g.posedit.pos
Interval(start=150550916, end=150550916, uncertain=False)

Start and end coordinates are polymorphic (can have multiple
representations). For genomic positions, these are instances of
:py:class:`SimplePosition`:

>>> var_g.posedit.pos.start
SimplePosition(base=150550916, uncertain=False)

For c. (cDNA) and r. (RNA) seqeunces, which have intron offsets and can be
measured from sequence start, CDS start, or CDS end (stop codon),
coordinates are more complex:

>>> var_c1.posedit.pos.start
BaseOffsetPosition(base=281, offset=0, datum=1, uncertain=False)

Either way, the sequence coordinate may be accessed via the ``base`` attribute:

>>> var_g.posedit.pos.start.base
150550916
>>> type(var_g.posedit.pos.start.base)
int


Map the genomic variant to another transcript
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

To map our genomic variant to another transcript, we need to provide a
transcript accession. One way to get those is to ask the data
provider:

>>> [ tx['ac'] for tx in hdp.get_tx_for_gene('MCL1') ]
['NM_021960.4', 'NM_182763.2', 'NM_001197320.1']

Let's map to the transcript for which this is an intronic variant.

>>> var_c2 = variantmapper.g_to_c(var_g,'NM_182763.2','GRCh37.p10')
>>> var_c2
SequenceVariant(ac=NM_182763.2, type=c, posedit=688+403C>T)
>>> var_c2.posedit.pos.start
BaseOffsetPosition(base=688, offset=403, datum=1, uncertain=False)

And, if we attempt to infer a protein consequence for this variant, we get
the expected uncertain interpretation:

>>> var_p2 = variantmapper.c_to_p(var_c2,None)
>>> var_p2
SequenceVariant(ac=NP_877495.1, type=p, posedit=?)
>>> str(var_p2)
'NP_877495.1:p.?'


