# -*- coding: utf-8 -*-
"""pyparsing-based grammar for HGVS variant strings.

Translates the OMeta grammar (formerly in hgvs/_data/hgvs.pymeta) into
pyparsing ParserElement objects.  Rules are built bottom-up and exposed
as attributes of the HGVSGrammar class.

cf. hgvs.pymeta for the original OMeta grammar and line references.
"""

import contextlib

import pyparsing as pp

import bioutils.sequences

import hgvs.edit
import hgvs.enums
import hgvs.hgvsposition
import hgvs.location
import hgvs.posedit
import hgvs.sequencevariant


@contextlib.contextmanager
def _hgvs_whitespace():
    """Build ParserElements with HGVS whitespace rules (i.e. none at all).

    CRITICAL: HGVS strings contain NO whitespace.

    pyparsing bakes DEFAULT_WHITE_CHARS into each element as it is constructed,
    so this must be in effect while the rules are built, not merely while they
    are used.

    The setting is process-global, and set_default_whitespace_chars() also
    rewrites pyparsing's builtin expressions in place, so it is restored on exit:
    leaving it applied would silently break unrelated pyparsing grammars in the
    same process.
    """
    prev = pp.ParserElement.DEFAULT_WHITE_CHARS
    pp.ParserElement.set_default_whitespace_chars("")
    try:
        yield
    finally:
        pp.ParserElement.set_default_whitespace_chars(prev)


class _NoneResult:
    """Sentinel for parse actions that must return Python None.

    pyparsing treats a parse-action return of None as "don't modify results",
    so we use this sentinel and unwrap it in HGVSGrammar.parse().
    """
    pass


_NONE_RESULT = _NoneResult()


def _s(val):
    """Unwrap a single-element ParseResults to a plain Python value."""
    if isinstance(val, pp.ParseResults):
        if len(val) == 1:
            return val[0]
        if len(val) == 0:
            return None
    return val


class HGVSGrammar:
    """All HGVS grammar rules as pyparsing ParserElements.

    Rules are built bottom-up:
      basic types -> locations -> edits -> posedits -> variants
    Each public attribute that is a ParserElement is collected into
    ``self.rules`` for name-based lookup.
    """

    def __init__(self):
        with _hgvs_whitespace():
            # Forward reference needed by dna_con/rna_con -> hgvs_position
            self.hgvs_position = pp.Forward()

            self._build_basic_types()
            self._build_locations()
            self._build_edits()
            self._build_posedits()
            self._build_typed_posedits()
            self._build_positions()
            self._build_variants()
            self._collect_rules()
            self._anchor_rules()

    def parse(self, rule_name, input_string):
        """Parse *input_string* using the named rule. Returns the domain object."""
        rule = self._anchored_rules[rule_name]
        result = rule.parse_string(input_string)
        rv = result[0]
        return None if isinstance(rv, _NoneResult) else rv

    # ------------------------------------------------------------------

    def _build_basic_types(self):
        """cf. hgvs.pymeta lines 164-218"""

        self.dna_iupac = pp.Char("ACGTRYMKWSBDHVNacgtrymkwsbdhvn")
        self.rna_iupac = pp.Char("ACGURYMKWSBDHVNacgurymkwsbdhvn")
        self.na_iupac = pp.Char("ACGTURYMKWSBDHVNacgturymkwsbdhvn")
        self.dna = self.dna_iupac
        self.rna = self.rna_iupac

        self.aa1 = pp.Char("ACDEFGHIKLMNPQRSTVWYBZXU")
        self.term1 = pp.Char("X*")
        self.pm = pp.Char("-+")

        # Not derived from bioutils.sequences.aa3_to_aa1_lut because that
        # table is missing Asx/Glx (ambiguous codes) and includes Ter
        # (which hgvs treats separately via term3).
        _aa3_list = [
            "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile",
            "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser",
            "Thr", "Val", "Trp", "Tyr", "Asx", "Glx", "Xaa", "Sec",
        ]
        self.aa3 = pp.MatchFirst([pp.Literal(aa) for aa in _aa3_list])
        self.term3 = pp.Literal("Ter")

        self.aa13 = self.aa3 | self.aa1
        self.aat1 = self.term1 | self.aa1
        self.aat3 = self.term3 | self.aa3
        self.aat13 = self.aat3 | self.aat1
        self.term13 = self.term3 | self.term1

        # --- Sequences ---
        self.dna_seq = pp.Combine(self.dna[1, ...])
        self.rna_seq = pp.Combine(self.rna[1, ...])
        self.aa1_seq = pp.Combine(self.aa1[1, ...])
        self.aa3_seq = pp.Combine(self.aa3[1, ...])
        self.aat1_seq = pp.Combine(self.term1) | pp.Combine(self.aa1[1, ...] + pp.Optional(self.term1))
        self.aat3_seq = pp.Combine(self.term3) | pp.Combine(self.aa3[1, ...] + pp.Optional(self.term3))
        self.aat13_seq = self.aat3_seq | self.aat1_seq
        self.aa13_seq = self.aa3_seq | self.aa1_seq

        # --- Numeric primitives ---
        _num_str = pp.Word(pp.nums)
        self.num = _num_str.copy().add_parse_action(lambda t: int(t[0]))
        self.nnum = pp.Combine(pp.Literal("-") + _num_str).add_parse_action(lambda t: int(t[0]))
        self.snum = pp.Combine(pp.Optional(self.pm) + _num_str).add_parse_action(lambda t: int(t[0]))
        self.base = self.snum.copy()
        self.offset = self.snum | pp.Empty().add_parse_action(pp.replace_with(0))

        # --- fs/ext helpers ---
        self.fsext_offset = (
            self.num
            | pp.Literal("?")
            | pp.Empty().add_parse_action(pp.replace_with(None))
        )

        self.aa13_fs = (
            pp.Suppress(self.term13) + self.fsext_offset
        ).add_parse_action(lambda t: t[0])

        self.fs = (
            pp.Suppress(pp.Literal("fs")) + (
                self.aa13_fs | pp.Empty().add_parse_action(pp.replace_with(None))
            )
        ).add_parse_action(lambda t: t[0])

        self.aa13_ext = (
            (self.term13 + self.fsext_offset).add_parse_action(
                lambda t: (t[0], t[1])
            )
            | ((self.aa13 | pp.Empty().add_parse_action(pp.replace_with(None))) + self.nnum).add_parse_action(
                lambda t: (t[0], t[1])
            )
        )

        self.ext = (
            pp.Suppress(pp.Literal("ext")) + (
                self.aa13_ext
                | pp.Empty().add_parse_action(lambda: [(None, None)])
            )
        ).add_parse_action(lambda t: t[0])

        # --- Accession ---
        self.accn = pp.Regex(r"[A-Za-z](?:[A-Za-z0-9]|[-_](?=[A-Za-z0-9]))*(?:\.[0-9]+)?")

        # --- Gene ---
        self.gene_symbol = pp.Regex(r"[A-Za-z](?:[A-Za-z0-9]|[-_](?=[A-Za-z0-9]))+")
        self.paren_gene = pp.Suppress("(") + self.gene_symbol + pp.Suppress(")")
        self.opt_gene_expr = self.paren_gene | pp.Empty().add_parse_action(pp.replace_with(None))

    def _build_locations(self):
        """cf. hgvs.pymeta lines 118-161"""

        # --- Definite positions ---
        self.def_c_pos = (
            (self.base + self.offset).add_parse_action(
                lambda t: hgvs.location.BaseOffsetPosition(t[0], t[1], datum=hgvs.enums.Datum.CDS_START)
            )
            | (pp.Suppress("*") + self.num + self.offset).add_parse_action(
                lambda t: hgvs.location.BaseOffsetPosition(t[0], t[1], datum=hgvs.enums.Datum.CDS_END)
            )
        )

        _def_gm_pos = (
            self.num | pp.Literal("?").add_parse_action(pp.replace_with(None))
        ).add_parse_action(
            lambda t: hgvs.location.SimplePosition(t[0])
        )
        self.def_g_pos = _def_gm_pos.copy()
        self.def_m_pos = _def_gm_pos.copy()

        self.def_n_pos = (self.base + self.offset).add_parse_action(
            lambda t: hgvs.location.BaseOffsetPosition(t[0], t[1], datum=hgvs.enums.Datum.SEQ_START)
        )

        self.def_p_pos = (
            (self.term13 | self.aa13) + self.num
        ).add_parse_action(
            lambda t: hgvs.location.AAPosition(t[1], bioutils.sequences.aa_to_aa1(t[0]))
        )

        self.def_r_pos = (self.base + self.offset).add_parse_action(
            lambda t: hgvs.location.BaseOffsetPosition(t[0], t[1], datum=hgvs.enums.Datum.SEQ_START)
        )

        # --- Positions (alias for definite) ---
        self.c_pos = self.def_c_pos
        self.g_pos = self.def_g_pos
        self.m_pos = self.def_m_pos
        self.n_pos = self.def_n_pos
        self.p_pos = self.def_p_pos
        self.r_pos = self.def_r_pos

        # --- Definite intervals ---
        self.def_g_interval = (
            (self.g_pos + pp.Suppress("_") + self.g_pos).add_parse_action(
                lambda t: hgvs.location.Interval(t[0], t[1])
            )
            | self.g_pos.copy().add_parse_action(
                lambda t: hgvs.location.Interval(t[0], None)
            )
        )

        self.def_m_interval = (
            (self.m_pos + pp.Suppress("_") + self.m_pos).add_parse_action(
                lambda t: hgvs.location.Interval(t[0], t[1])
            )
            | self.m_pos.copy().add_parse_action(
                lambda t: hgvs.location.Interval(t[0], None)
            )
        )

        self.def_p_interval = (
            (self.p_pos + pp.Suppress("_") + self.p_pos).add_parse_action(
                lambda t: hgvs.location.Interval(t[0], t[1])
            )
            | self.p_pos.copy().add_parse_action(
                lambda t: hgvs.location.Interval(t[0], None)
            )
        )

        self.def_r_interval = (
            (self.r_pos + pp.Suppress("_") + self.r_pos).add_parse_action(
                lambda t: hgvs.location.Interval(t[0], t[1])
            )
            | self.r_pos.copy().add_parse_action(
                lambda t: hgvs.location.Interval(t[0], None)
            )
        )

        self.def_c_interval = (
            (self.c_pos + pp.Suppress("_") + self.c_pos).add_parse_action(
                lambda t: hgvs.location.BaseOffsetInterval(t[0], t[1])
            )
            | self.c_pos.copy().add_parse_action(
                lambda t: hgvs.location.BaseOffsetInterval(t[0], None)
            )
        )

        self.def_n_interval = (
            (self.n_pos + pp.Suppress("_") + self.n_pos).add_parse_action(
                lambda t: hgvs.location.BaseOffsetInterval(t[0], t[1])
            )
            | self.n_pos.copy().add_parse_action(
                lambda t: hgvs.location.BaseOffsetInterval(t[0], None)
            )
        )

        # --- Uncertain genomic intervals ---
        self.uncertain_g_interval = (
            (pp.Suppress("(") + self.def_g_interval + pp.Suppress(")") +
             pp.Suppress("_") +
             pp.Suppress("(") + self.def_g_interval + pp.Suppress(")")).add_parse_action(
                lambda t: hgvs.location.Interval(start=t[0]._set_uncertain(), end=t[1]._set_uncertain())
            )
            | (self.def_g_interval + pp.Suppress("_") +
               pp.Suppress("(") + self.def_g_interval + pp.Suppress(")")).add_parse_action(
                lambda t: hgvs.location.Interval(start=t[0], end=t[1]._set_uncertain())
            )
            | (pp.Suppress("(") + self.def_g_interval + pp.Suppress(")") +
               pp.Suppress("_") + self.def_g_interval).add_parse_action(
                lambda t: hgvs.location.Interval(start=t[0]._set_uncertain(), end=t[1])
            )
            | (pp.Suppress("(") + self.def_g_interval + pp.Suppress(")")).add_parse_action(
                lambda t: t[0]._set_uncertain()
            )
        )

        # --- Potentially uncertain intervals ---
        self.c_interval = (
            self.def_c_interval
            | (pp.Suppress("(") + self.def_c_interval + pp.Suppress(")")).add_parse_action(
                lambda t: t[0]._set_uncertain()
            )
        )
        self.g_interval = self.uncertain_g_interval | self.def_g_interval
        self.m_interval = (
            self.def_m_interval
            | (pp.Suppress("(") + self.def_m_interval + pp.Suppress(")")).add_parse_action(
                lambda t: t[0]._set_uncertain()
            )
        )
        self.n_interval = (
            self.def_n_interval
            | (pp.Suppress("(") + self.def_n_interval + pp.Suppress(")")).add_parse_action(
                lambda t: t[0]._set_uncertain()
            )
        )
        self.p_interval = (
            self.def_p_interval
            | (pp.Suppress("(") + self.def_p_interval + pp.Suppress(")")).add_parse_action(
                lambda t: t[0]._set_uncertain()
            )
        )
        self.r_interval = (
            self.def_r_interval
            | (pp.Suppress("(") + self.def_r_interval + pp.Suppress(")")).add_parse_action(
                lambda t: t[0]._set_uncertain()
            )
        )

    def _build_edits(self):
        """cf. hgvs.pymeta lines 76-116"""

        _num_str = pp.Word(pp.nums)

        # --- DNA edits ---
        self.dna_ident = (
            pp.Combine(self.dna[...]) + pp.Suppress(pp.Literal("="))
        ).add_parse_action(
            lambda t: hgvs.edit.NARefAlt(ref=t[0], alt=t[0])
        )

        self.dna_subst = (
            self.dna + pp.Suppress(pp.Literal(">")) + self.dna
        ).add_parse_action(
            lambda t: hgvs.edit.NARefAlt(ref=t[0], alt=t[1])
        )

        self.dna_delins = (
            pp.Suppress(pp.Literal("del")) +
            pp.Combine(_num_str | self.dna[...]) +
            pp.Suppress(pp.Literal("ins")) +
            pp.Combine(self.dna[1, ...])
        ).add_parse_action(
            lambda t: hgvs.edit.NARefAlt(ref=t[0], alt=t[1])
        )

        self.dna_del = (
            pp.Suppress(pp.Literal("del")) +
            pp.Combine(_num_str | self.dna[...])
        ).add_parse_action(
            lambda t: hgvs.edit.NARefAlt(ref=t[0], alt=None)
        )

        self.dna_ins = (
            pp.Suppress(pp.Literal("ins")) + pp.Combine(self.dna[1, ...])
        ).add_parse_action(
            lambda t: hgvs.edit.NARefAlt(ref=None, alt=t[0])
        )

        self.dna_dup = (
            pp.Suppress(pp.Literal("dup")) + pp.Combine(self.dna[...])
        ).add_parse_action(
            lambda t: hgvs.edit.Dup(ref=t[0])
        )

        self.dna_inv = (
            pp.Suppress(pp.Literal("inv")) + pp.Combine(_num_str | self.dna[...])
        ).add_parse_action(
            lambda t: hgvs.edit.Inv(ref=None)
        )

        self.dna_con = (
            pp.Suppress(pp.Literal("con")) + self.hgvs_position
        ).add_parse_action(
            lambda t: hgvs.edit.Conv(from_ac=t[0].ac, from_type=t[0].type, from_pos=t[0].pos)
        )

        self.dna_copy = (
            pp.Suppress(pp.Literal("copy")) + self.num
        ).add_parse_action(
            lambda t: hgvs.edit.NACopy(copy=t[0])
        )

        self.dna_edit = (
            self.dna_ident | self.dna_subst | self.dna_delins | self.dna_ins
            | self.dna_del | self.dna_dup | self.dna_inv | self.dna_con | self.dna_copy
        )

        self.dna_edit_mu = (
            self.dna_edit
            | (pp.Suppress("(") + self.dna_edit + pp.Suppress(")")).add_parse_action(
                lambda t: t[0]._set_uncertain()
            )
        )

        # --- RNA edits ---
        self.rna_ident = (
            pp.Combine(self.rna[...]) + pp.Suppress(pp.Literal("="))
        ).add_parse_action(
            lambda t: hgvs.edit.NARefAlt(ref=t[0], alt=t[0])
        )

        self.rna_subst = (
            self.rna + pp.Suppress(pp.Literal(">")) + self.rna
        ).add_parse_action(
            lambda t: hgvs.edit.NARefAlt(ref=t[0], alt=t[1])
        )

        self.rna_delins = (
            pp.Suppress(pp.Literal("del")) +
            pp.Combine(_num_str | self.rna[...]) +
            pp.Suppress(pp.Literal("ins")) +
            pp.Combine(self.rna[1, ...])
        ).add_parse_action(
            lambda t: hgvs.edit.NARefAlt(ref=t[0], alt=t[1])
        )

        self.rna_del = (
            pp.Suppress(pp.Literal("del")) +
            pp.Combine(_num_str | self.rna[...])
        ).add_parse_action(
            lambda t: hgvs.edit.NARefAlt(ref=t[0], alt=None)
        )

        self.rna_ins = (
            pp.Suppress(pp.Literal("ins")) + pp.Combine(self.rna[1, ...])
        ).add_parse_action(
            lambda t: hgvs.edit.NARefAlt(ref=None, alt=t[0])
        )

        self.rna_dup = (
            pp.Suppress(pp.Literal("dup")) + pp.Combine(self.rna[...])
        ).add_parse_action(
            lambda t: hgvs.edit.Dup(ref=t[0])
        )

        self.rna_inv = (
            pp.Suppress(pp.Literal("inv")) + pp.Combine(_num_str | self.rna[...])
        ).add_parse_action(
            lambda t: hgvs.edit.Inv(ref=None)
        )

        self.rna_con = (
            pp.Suppress(pp.Literal("con")) + self.hgvs_position
        ).add_parse_action(
            lambda t: hgvs.edit.Conv(from_ac=t[0].ac, from_type=t[0].type, from_pos=t[0].pos)
        )

        self.rna_edit = (
            self.rna_ident | self.rna_subst | self.rna_delins | self.rna_ins
            | self.rna_del | self.rna_dup | self.rna_inv | self.rna_con
        )

        self.rna_edit_mu = (
            self.rna_edit
            | (pp.Suppress("(") + self.rna_edit + pp.Suppress(")")).add_parse_action(
                lambda t: t[0]._set_uncertain()
            )
        )

        # --- Protein edits ---
        self.pro_subst = (
            self.aat13 | pp.Literal("?")
        ).add_parse_action(
            lambda t: hgvs.edit.AASub(ref="", alt=t[0])
        )

        self.pro_delins = (
            pp.Suppress(pp.Literal("delins")) + self.aat13_seq
        ).add_parse_action(
            lambda t: hgvs.edit.AARefAlt(ref="", alt=t[0])
        )

        self.pro_del = pp.Literal("del").add_parse_action(
            lambda: hgvs.edit.AARefAlt(ref="", alt=None)
        )

        self.pro_ins = (
            pp.Suppress(pp.Literal("ins")) + self.aat13_seq
        ).add_parse_action(
            lambda t: hgvs.edit.AARefAlt(ref=None, alt=t[0])
        )

        self.pro_dup = pp.Literal("dup").add_parse_action(
            lambda: hgvs.edit.Dup(ref="")
        )

        self.pro_fs = (
            (self.aat13 | pp.Empty().add_parse_action(pp.replace_with("")))
            + self.fs
        ).add_parse_action(
            lambda t: hgvs.edit.AAFs(ref="", alt=t[0], length=t[1])
        )

        self.pro_ext = (
            pp.Optional(self.aat13, default=None) + self.ext
        ).add_parse_action(
            lambda t: hgvs.edit.AAExt(
                ref="", alt=t[0], aaterm=t[1][0], length=t[1][1]
            )
        )

        self.pro_ident = pp.Literal("=").add_parse_action(
            lambda: hgvs.edit.AARefAlt(ref="", alt="")
        )

        self.pro_edit = (
            self.pro_fs | self.pro_ext | self.pro_subst | self.pro_delins
            | self.pro_ins | self.pro_del | self.pro_dup | self.pro_ident
        )

        self.pro_edit_mu = (
            self.pro_edit
            | (pp.Suppress("(") + self.pro_edit + pp.Suppress(")")).add_parse_action(
                lambda t: t[0]._set_uncertain()
            )
        )

    def _build_posedits(self):
        """cf. hgvs.pymeta lines 56-73"""

        self.c_posedit = (self.c_interval + self.dna_edit).add_parse_action(
            lambda t: hgvs.posedit.PosEdit(pos=t[0], edit=t[1])
        )
        self.g_posedit = (self.g_interval + self.dna_edit).add_parse_action(
            lambda t: hgvs.posedit.PosEdit(pos=t[0], edit=t[1])
        )
        self.m_posedit = (self.m_interval + self.dna_edit).add_parse_action(
            lambda t: hgvs.posedit.PosEdit(pos=t[0], edit=t[1])
        )
        self.n_posedit = (self.n_interval + self.dna_edit).add_parse_action(
            lambda t: hgvs.posedit.PosEdit(pos=t[0], edit=t[1])
        )

        self.r_posedit = (
            (self.r_interval + self.rna_edit).add_parse_action(
                lambda t: hgvs.posedit.PosEdit(pos=t[0], edit=t[1])
            )
            | (pp.Suppress("(") + self.r_interval + self.rna_edit + pp.Suppress(")")).add_parse_action(
                lambda t: hgvs.posedit.PosEdit(pos=t[0], edit=t[1], uncertain=True)
            )
        )

        self.p_posedit_special = (
            pp.Literal("=").add_parse_action(
                lambda: hgvs.posedit.PosEdit(pos=None, edit="=", uncertain=False)
            )
            | (pp.Suppress("(") + pp.Suppress("=") + pp.Suppress(")")).add_parse_action(
                lambda: hgvs.posedit.PosEdit(pos=None, edit="=", uncertain=True)
            )
            | (pp.Literal("0") + pp.Suppress("?")).add_parse_action(
                lambda: hgvs.posedit.PosEdit(pos=None, edit="0", uncertain=True)
            )
            | pp.Literal("0").add_parse_action(
                lambda: hgvs.posedit.PosEdit(pos=None, edit="0", uncertain=False)
            )
            | pp.Literal("?").add_parse_action(lambda: _NONE_RESULT)
        )

        self.p_posedit = (
            (self.p_interval + self.pro_edit).add_parse_action(
                lambda t: hgvs.posedit.PosEdit(pos=t[0], edit=t[1])
            )
            | (pp.Suppress("(") + self.p_interval + self.pro_edit + pp.Suppress(")")).add_parse_action(
                lambda t: hgvs.posedit.PosEdit(pos=t[0], edit=t[1], uncertain=True)
            )
            | self.p_posedit_special
        )

    def _build_typed_posedits(self):
        """cf. hgvs.pymeta lines 42-52"""

        def _make_typed(type_char, posedit_rule):
            return (
                pp.Literal(type_char) + pp.Suppress(".") + posedit_rule
            ).add_parse_action(
                lambda t: hgvs.sequencevariant.SequenceVariant(
                    ac=None, type=t[0], posedit=None if isinstance(t[1], _NoneResult) else t[1]
                )
            )

        self.c_typed_posedit = _make_typed("c", self.c_posedit)
        self.g_typed_posedit = _make_typed("g", self.g_posedit)
        self.m_typed_posedit = _make_typed("m", self.m_posedit)
        self.n_typed_posedit = _make_typed("n", self.n_posedit)
        self.p_typed_posedit = _make_typed("p", self.p_posedit)
        self.r_typed_posedit = _make_typed("r", self.r_posedit)

    def _build_positions(self):
        """cf. hgvs.pymeta lines 28-38"""

        def _make_hgvs_pos(type_char, interval_rule):
            return (
                self.accn + self.opt_gene_expr + pp.Suppress(":") +
                pp.Literal(type_char) + pp.Suppress(".") + interval_rule
            ).add_parse_action(
                lambda t: hgvs.hgvsposition.HGVSPosition(ac=t[0], gene=_s(t[1]), type=t[2], pos=t[3])
            )

        self.c_hgvs_position = _make_hgvs_pos("c", self.c_interval)
        self.g_hgvs_position = _make_hgvs_pos("g", self.g_interval)
        self.m_hgvs_position = _make_hgvs_pos("m", self.m_interval)
        self.n_hgvs_position = _make_hgvs_pos("n", self.n_interval)
        self.p_hgvs_position = _make_hgvs_pos("p", self.p_interval)
        self.r_hgvs_position = _make_hgvs_pos("r", self.r_interval)

        # Resolve the Forward reference
        self.hgvs_position <<= (
            self.g_hgvs_position | self.m_hgvs_position | self.c_hgvs_position
            | self.n_hgvs_position | self.r_hgvs_position | self.p_hgvs_position
        )

    def _build_variants(self):
        """cf. hgvs.pymeta lines 14-24"""

        def _make_variant(type_char, posedit_rule):
            return (
                self.accn + self.opt_gene_expr + pp.Suppress(":") +
                pp.Literal(type_char) + pp.Suppress(".") + posedit_rule
            ).add_parse_action(
                lambda t: hgvs.sequencevariant.SequenceVariant(
                    ac=t[0], gene=_s(t[1]), type=t[2],
                    posedit=None if isinstance(t[3], _NoneResult) else t[3]
                )
            )

        self.c_variant = _make_variant("c", self.c_posedit)
        self.g_variant = _make_variant("g", self.g_posedit)
        self.m_variant = _make_variant("m", self.m_posedit)
        self.n_variant = _make_variant("n", self.n_posedit)
        self.p_variant = _make_variant("p", self.p_posedit)
        self.r_variant = _make_variant("r", self.r_posedit)

        self.hgvs_variant = (
            self.g_variant | self.m_variant | self.c_variant
            | self.n_variant | self.r_variant | self.p_variant
        )

    def _collect_rules(self):
        """Build name->rule dict from all instance attributes that are ParserElements."""
        self.rules = {}
        for name in dir(self):
            if name.startswith("_"):
                continue
            obj = getattr(self, name)
            if isinstance(obj, pp.ParserElement):
                self.rules[name] = obj

    def _anchor_rules(self):
        """Pair each rule with an end-of-string anchor for whole-input matching.

        Must be called while the no-whitespace default is in effect (see
        _hgvs_whitespace). parse_string(parse_all=True) is not used because it
        builds its own Empty() + StringEnd() at *parse* time, which picks up
        whatever whitespace default is globally in effect then and would let
        trailing whitespace through.
        """
        self._anchored_rules = {
            name: rule + pp.StringEnd() for name, rule in self.rules.items()
        }


# <LICENSE>
# Copyright 2018 HGVS Contributors (https://github.com/biocommons/hgvs)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
