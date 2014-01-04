Genome, transcript, and protein sequence variants are typically reported
using the `mutation nomenclature ("mutnomen") recommendations
<http://www.hgvs.org/mutnomen/>`_ provided by the `Human Genome Variation
Society (HGVS) <http://www.hgvs.org/>`_. This standard was formulated in
an era of specialized sequencing and cytogenetic analyses, long before
high-throughput sequencing annotation was envisioned.  Unfortunately, the
complexity of biological phenomena and the breadth of the standard makes
it difficult to implement the standard in software.

This package, ``hgvs``, is an easy-to-use Python library for parsing,
representing, formatting, and mapping variants between genome, transcript,
and protein sequences.  The current implementation handles most (but not
all) of the standard for precisely defined sequence variants.  The intent
is to centralize the subset of HGVS variant manipulation that is routinely
used in modern, high-throughput sequencing analysis.

