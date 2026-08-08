.. _shift_over_boundary:

Boundary-Spanning Variants: ``shift_over_boundary``
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Background
@@@@@@@@@@

The HGVS 3' shifting rule is a nomenclature convention — all variants are shifted as far as
possible toward the 3' end of the transcript. This is arbitrary from a biological standpoint:
biology does not follow the 3' rule, and in some cases strict application of it causes
information to be lost.

Consider a positive-strand gene where the first 8 bases of an intron are duplicated::

    c.100+1_100+8dup   (intronic positions, both have an offset)

Under standard mapping this yields ``p.?`` because both positions are intronic and no protein
consequence can be determined. But consider the biology: after the duplication there are now two
potential splice sites on the left side of the intron. The most parsimonious assumption is that
splicing uses the *inner* (original) splice site, leaving the extra material inside the coding
sequence — which would produce a frameshift. The entire intronic sequence downstream of the
duplication remains intact, so the splice site is preserved.

This feature was motivated in part by **FLT3 Internal Tandem Duplications (ITDs)**, a variant
class in acute myeloid leukemia. A fraction of FLT3 ITDs span the exon/intron
boundary and were annotated as unclear protein effect (returning ``p.?``) instead of producing a specific
protein consequence.


How it works
@@@@@@@@@@@@

When ``shift_over_boundary`` is enabled and a ``c_to_p`` mapping returns ``p.?`` for an ``ins``
or ``dup`` variant at an intron/exon boundary, the mapper tries shifting the variant in the
*reverse* direction (i.e., opposite to the HGVS 3' rule). If a shifted representation exists that
places the variant entirely within an exon or entirely within an intron, and the splice site and
splice region remain intact, that shifted form is used for protein consequence prediction.

The HGVS notation of the input ``c.`` variant is **not changed** — only the internal
representation used for ``p.`` prediction is affected. This means the output remains
nomenclature-compliant.


Configuration
@@@@@@@@@@@@@

Two settings control this behaviour. They can be set per-mapper instance in code, or globally
via ``hgvs.global_config``.


``shift_over_boundary`` (bool, default: ``False``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Master switch. Disabled by default to preserve backward compatibility with existing code and
tests.

.. code-block:: python

    am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37")
    am.shift_over_boundary = True


``shift_over_boundary_preference`` (enum, default: ``DEFAULT``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When a variant can be shifted to either side of the boundary, this controls which side is
preferred.

+-------------+------------------------------------------------------------------------------------+
| Value       | Behaviour                                                                          |
+=============+====================================================================================+
| ``DEFAULT`` | Use the representation that results from standard HGVS rules; only use a shifted  |
|             | form if it uniquely resolves the variant                                           |
+-------------+------------------------------------------------------------------------------------+
| ``EXON``    | Prefer the exonic representation — yields a specific protein consequence when the  |
|             | splice region is intact                                                            |
+-------------+------------------------------------------------------------------------------------+
| ``INTRON``  | Prefer the intronic representation — yields ``p.?`` (protein consequence unknown) |
+-------------+------------------------------------------------------------------------------------+

.. code-block:: python

    from hgvs.enums import ShiftOverBoundaryPreference

    am.shift_over_boundary = True
    am.shift_over_boundary_preference = ShiftOverBoundaryPreference.EXON


Example: 2-base dup at an intron/exon boundary (*EZH2*)
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

``NM_004456.4:c.2196-1_2196dup`` duplicates the last intronic base and first exonic base of
*EZH2*. It is a compact, clear illustration of the feature.

.. code-block:: python

    import hgvs.assemblymapper
    import hgvs.parser
    from hgvs.enums import ShiftOverBoundaryPreference

    hp = hgvs.parser.Parser()
    am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37")

    var_g = hp.parse("NC_000007.13:g.148504799_148504800dup")

    # Default behaviour (shift_over_boundary disabled)
    am.shift_over_boundary = False
    print(am.g_to_c(var_g, "NM_004456.4"))   # NM_004456.4:c.2196-2_2196-1dup
    print(am.c_to_p(am.g_to_c(var_g, "NM_004456.4")))
    # NP_004447.2:p.?

    # Shift enabled, prefer exonic representation → specific protein consequence
    am.shift_over_boundary = True
    am.shift_over_boundary_preference = ShiftOverBoundaryPreference.EXON
    print(am.c_to_p(am.g_to_c(var_g, "NM_004456.4")))
    # NP_004447.2:p.(Tyr733AspfsTer8)

    # Shift enabled, prefer intronic representation → unknown consequence
    am.shift_over_boundary_preference = ShiftOverBoundaryPreference.INTRON
    print(am.c_to_p(am.g_to_c(var_g, "NM_004456.4")))
    # NP_004447.2:p.?


Pretty-print view
~~~~~~~~~~~~~~~~~

Using ``PrettyPrint`` makes the effect visible in genomic context. The ``hgvs_c`` and ``hgvs_g``
representations are identical in both cases — only ``hgvs_p`` changes.

**Default (** ``shift_over_boundary=False`` **) — protein consequence unknown:**

.. code-block:: text

    hgvs_g    : NC_000007.13:g.148504799_148504800dup
    hgvs_c    : NM_004456.4:c.2196-2_2196-1dup
    hgvs_p    : NP_004447.2:p.?
              : 148,504,810         148,504,790
    chrom pos : |    .    |    .    |
    seq    <- : GACAACAAAGTCTATGTCGGTCC
    seq    -> : CTGTTGTTTCAGATACAGCCAGG
    region    :           |-|
    tx seq -> :             ATACAGCCAGG
    tx pos    :                 |    .
              :                 2200
    aa seq -> :             gTyrSerGlnA

**EXON preference (** ``shift_over_boundary=True`` **,** ``ShiftOverBoundaryPreference.EXON`` **) — frameshift predicted:**

.. code-block:: text

    hgvs_g    : NC_000007.13:g.148504799_148504800dup
    hgvs_c    : NM_004456.4:c.2196-2_2196-1dup
    hgvs_p    : NP_004447.2:p.(Tyr733AspfsTer8)
              : 148,504,810         148,504,790
    chrom pos : |    .    |    .    |
    seq    <- : GACAACAAAGTCTATGTCGGTCC
    seq    -> : CTGTTGTTTCAGATACAGCCAGG
    region    :           |-|
    tx seq -> :             ATACAGCCAGG
    tx pos    :                 |    .
              :                 2200
    aa seq -> :             gTyrSerGlnA

The ``hgvs_p`` line changes from ``p.?`` to ``p.(Tyr733AspfsTer8)`` — a frameshift at Tyr733
with a stop codon 8 residues later.


Example: FLT3 ITD spanning an intron/exon boundary (longer insertion)
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

For a longer insertion, see a FLT3 ITD case:

.. code-block:: python

    var_c = hp.parse("NM_004119.2:c.1837+21_1837+22insCGAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAA")

    am.shift_over_boundary = True
    am.shift_over_boundary_preference = ShiftOverBoundaryPreference.EXON
    print(am.c_to_p(var_c))
    # NP_004110.2:p.(Lys614_Val615insAsnGlyMetCysGlnThrArgGluTyrGluTyrAspLeuLysTrpGluPheProArgGluAsnLeuGluPheGlyLys)


When does shifting apply?
@@@@@@@@@@@@@@@@@@@@@@@@@

The mapper only attempts to shift when **all** of the following are true:

1. ``shift_over_boundary = True``
2. The variant type is ``ins`` or ``dup``
3. The variant region is intronic or exonic (the variant is at the boundary, not disrupting core coding sequence)
4. A shifted representation exists that moves the variant to the other side of the boundary
5. The splice site and splice region remain intact after the shift

If no valid shift is possible — for example, the variant genuinely disrupts the splice site, or
the shifted form produces an invalid interval — the original position is used and ``p.?`` is
returned regardless of preference.


Relationship to ``ins_at_boundary_is_intronic``
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

The existing ``ins_at_boundary_is_intronic`` setting governs how a single-position insertion
*exactly at* the exon/intron junction is classified. ``shift_over_boundary`` is complementary:
it applies to variants that already have an unambiguous HGVS position but whose duplicated or
inserted sequence happens to span the boundary, making the protein consequence undetermined under
strict HGVS rules.


Clinical significance
@@@@@@@@@@@@@@@@@@@@@

This feature was motivated by *FLT3* Internal Tandem Duplications (ITDs) in acute myeloid
leukemia (AML). FLT3-ITDs vary in length and position; boundary-spanning cases have been
documented in sequencing studies. Without ``shift_over_boundary``, these variants return ``p.?``,
which could cause them to be silently dropped in downstream annotation pipelines.


See also
@@@@@@@@

- GitHub issue `#714 <https://github.com/biocommons/hgvs/issues/714>`_ — original bug report
- GitHub PR `#719 <https://github.com/biocommons/hgvs/pull/719>`_ — implementation by Brendan O'Donnell (genomoncology)
- :class:`hgvs.enums.ShiftOverBoundaryPreference`
- :class:`hgvs.assemblymapper.AssemblyMapper`
- :class:`hgvs.variantmapper.VariantMapper`
