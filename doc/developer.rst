=====================
Developer Information
=====================

The Human Genome Variation Society, `HGVS`_, produced the `HGVS
Recommendations`_ for the representation of genome, transcript, and
protein, variation.  The following is a *brief* overview of the the
standard and especially the subset parsed by this package.

A simple (not compound) variant is composed of the following general
elements:

  <seqref> : <type> . <position> <edit>

  Example:
  NC_000006.1:g.13403110A>T
  NM_004006.1:c.3G>T
  NR_001234.5:r.15delA
  NP_009876.5:p.Trp26Cys

The present parser requires all components. Gene names are not accepted,
either in lieu of or in addition to the sequence name (the Recommendations
permit both). The recommendation uses parenthesis to express uncertainty
in position, edit, and other details; this is not supported.


Classes
-------


