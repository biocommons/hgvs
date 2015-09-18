"""Internal utility to display textual representation of variant context

Something like this:

(default-2.7)snafu$ ./sbin/hgvs-shell 
In [1]: import hgvs.utils.context as huc; h = 'NM_000399.3:n.2689_2690insA'; v = hp.parse_hgvs_variant(h); vn = hn.normalize(v)
In [2]: reload(huc); print(huc.variant_context_w_alignment(evm, vn))
                                              v                               NC_000010.10:g.64572045dupT
NC_000010.10 g 64572025 > ACTCAGGGAGTGATTTTTTTTCTCCATAATAAGGCAACCCA          > 64572065 NC_000010.10:g.64572045dupT
NC_000010.10 g 64572025 < TGAGTCCCTCACTAAAAAAAAGAGGTATTATTCCGTTGGGT          < 64572065 NC_000010.10:g.64572045dupT
                          |||||||||||||-|||||||||||||||||||||||||||          13=1D27= 
NM_000399.3  n     2670 < TGAGTCCCTCACT-AAAAAAAGAGGTATTATTCCGTTGGGT          <     2709 NM_000399.3:n.2696dupA
NM_000399.3  c      902 <                                                    <      941 NM_000399.3:c.*928dupA

"""

import re

from bioutils.sequences import complement, reverse_complement
from hgvs.location import Interval, SimplePosition

def full_house(evm, var, tx_ac=None):
    if var.type == 'g':
        var_g = var
        if tx_ac is None:
            rtx = evm.relevant_transcripts(var)
            if len(rtx) == 0:
                raise RuntimeError("no relevant transcripts for {var.ac}".format(var=var))
            if len(rtx) > 1:
                raise RuntimeError("{n} relevant transcripts for {var.ac}; you need to pick one".format(n=len(rtx),var=var))
            tx_ac = rtx[0]
        var_n = evm.g_to_n(var_g, tx_ac)
        var_c = evm.n_to_c(var_n)
    elif var.type == 'n':
        var_n = var
        var_g = evm.n_to_g(var_n)
        var_c = evm.n_to_c(var_n)
    elif var.type == 'c':
        var_c = var
        var_g = evm.c_to_g(var_c)
        var_n = evm.c_to_n(var_c)
    var_p = evm.c_to_p(var_c)
    return {'g': var_g, 'c': var_c, 'n': var_n, 'p': var_p}


def variant_context(evm, var, margin=20):
    span = _span(var, margin)
    span_g = _ival_to_span(fh['g'])
    span_g[0] -= margin
    span_g[1] += margin
    return '\n'.join([
        seq_line_fmt(var=var, span=span, content=evm.hdp.fetch_seq(var.ac, *span), post=''),
        pointer_line(var, span)
        ])


def variant_context_w_alignment(evm, var, margin=20, tx_ac=None):
    """This module is experimental. It requires the uta_align package from pypi."""
    from uta_align.align.algorithms import align, cigar_alignment

    fh = full_house(evm, var, tx_ac=tx_ac)
    tm = evm._fetch_TranscriptMapper(fh['n'].ac, fh['g'].ac, evm.alt_aln_method)
    strand = tm.strand
    span_g = _ival_to_span(fh['g'].posedit.pos)
    span_g = (span_g[0]-margin, span_g[1]+margin)
    ival_g = Interval(SimplePosition(span_g[0]), SimplePosition(span_g[1]))
    ival_n = tm.g_to_n(ival_g)
    assert ival_n.start.offset == 0 and ival_n.end.offset == 0, "limited to coding variants"
    span_n = _ival_to_span(ival_n)
    ival_c = tm.g_to_c(ival_g)
    span_c = _ival_to_span(ival_c)
    seq_gt = evm.hdp.fetch_seq(fh['g'].ac, span_g[0]-1, span_g[1])
    seq_gb = complement(seq_gt)
    seq_n = evm.hdp.fetch_seq(fh['n'].ac, span_n[0]-1, span_n[1])
    if strand == 1:
        a = align(seq_gt, seq_n, b'global', extended_cigar=True)
    else:
        seq_n = ''.join(reversed(seq_n))
        a = align(seq_gb, seq_n, b'global', extended_cigar=True)

    aseq_gt, _      = cigar_alignment(seq_gt, a.query, a.cigar, hide_match=False)
    aseq_gb, aseq_n = cigar_alignment(seq_gb, a.query, a.cigar, hide_match=False)
    aln_str = _reformat_aln_str(cigar_alignment(a.ref, a.query, a.cigar, hide_match=True)[1])
    s_dir = '>' if strand == 1 else '<'

    lines = [
        [1, 0, seq_line_fmt(var=fh['c'], span=span_c, content='',      dir=s_dir),              ],
        [2, 0, seq_line_fmt(var=fh['n'], span=span_n, content=aseq_n,  dir=s_dir),              ],
        [3, 0, _line_fmt.format(pre='', content=aln_str, post=a.cigar.to_string(), comment=''), ],
        [4, 1, seq_line_fmt(var=fh['g'], span=span_g, content=aseq_gt, dir='>'),                ],
        [4, 2, seq_line_fmt(var=fh['g'], span=span_g, content=aseq_gb, dir='<'),                ],
        [5, 0, pointer_line(var=fh['g'], span=span_g),                                          ],
        ]
    if strand == -1:
        lines.sort(key = lambda e: (-e[0], e[1]))
    return '\n'.join(r[2] for r in lines)


def _ival_to_span(ival):
    return (ival.start.base, ival.end.base)


def _reformat_aln_str(aln_str):
    return re.sub('[ACGT]', ' ', aln_str.replace('.', '|'))

# pre=[ac c s] d content d post=[end] comment
_line_fmt = "{pre:>30s} {content:45s} {post} {comment}"
_pre_fmt = "{ac:12s} {type:1s} {s:10d} {dir:1s}"
_post_fmt = "{dir:1s} {e:8d}"

def seq_line_fmt(var, span, content, dir=''):
    return _line_fmt.format(
        pre=_pre_fmt.format(ac=var.ac, type=var.type, s=span[0], dir=dir),
        content=content,
        post=_post_fmt.format(dir=dir, e=span[1]),
        comment=str(var))
def pointer_line(var, span):
    s0 = span[0]
    o = var.posedit.pos.start.base - s0
    l = var.posedit.pos.end.base - var.posedit.pos.start.base + 1
    if var.posedit.edit.type == 'ins':
        p = ' '*o + '><'
    else:
        p = ' '*o + '*'*l
    return _line_fmt.format(pre='', content=p, post='', comment=str(var))

