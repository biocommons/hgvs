# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

#
# Utility class for creating an hgvsp SequenceVariant object,
# given a transcript with variants applied.
# Used in hgvsc to hgvsp conversion.
#
import collections
import difflib

import hgvs
import hgvs.edit
import hgvs.exceptions
import hgvs.location
import hgvs.posedit
import hgvs.utils
import hgvs.variant

DBG = False


class AltSeqToHgvsp(object):
    def __init__(self, ref_data, alt_data):
        """Constructor

        :param ref_data: reference transcript record
        :type ref_data: recordtype
        :param alt_data: alt transcript record
        :type ref_data: recordtype
        """
        self._ref_data = ref_data
        self._alt_data = alt_data
        self._protein_accession = self._ref_data.protein_accession
        self._ref_seq = self._ref_data.aa_sequence
        self._alt_seq = self._alt_data.aa_sequence
        self._is_frameshift = self._alt_data.is_frameshift
        self._frameshift_start = self._alt_data.frameshift_start
        self._is_substitution = self._alt_data.is_substitution
        self._is_ambiguous = self._alt_data.is_ambiguous

        if DBG:
            print("len ref seq:{} len alt seq:{}".format(len(self._ref_seq), len(self._alt_seq)))
            print("fs start:{} protein ac:{}".format(self._frameshift_start, self._protein_accession))
            print(self._ref_seq)
            print(self._alt_seq)
            print("aa variant start: {}".format(self._alt_data.variant_start_aa))
            print(self._ref_data.transcript_sequence)
            print(self._alt_data.transcript_sequence)

    def build_hgvsp(self):
        """Compare two amino acid sequences; generate an hgvs tag from the output

        :return list of variants in sequence order
        :rtype list of dict
        """

        variants = []

        if not self._is_ambiguous and len(self._alt_seq) > 0:

            do_delins = True
            if self._ref_seq == self._alt_seq:
                do_delins = False
            elif self._is_substitution:
                if len(self._ref_seq) == len(self._alt_seq):
                    diff_pos = [(i, self._ref_seq[i], self._alt_seq[i]) for i in xrange(len(self._ref_seq))
                                if self._ref_seq[i] != self._alt_seq[i]]
                    if len(diff_pos) == 1:
                        (start, deletion, insertion) = diff_pos[0]
                        variants.append({"start": start + 1, "ins": insertion, "del": deletion})
                        do_delins = False

                elif self._alt_seq[self._alt_data.variant_start_aa - 1] == "*" and \
                                self._ref_seq[self._alt_data.variant_start_aa - 1] != "*":    # introduced stop codon
                    deletion = self._ref_seq[self._alt_data.variant_start_aa - 1:]
                    variants.append({"start": self._alt_data.variant_start_aa, "ins": "*", "del": deletion})
                    do_delins = False

            if do_delins:
                if self._alt_data.is_frameshift:
                    start = self._alt_data.variant_start_aa - 1
                    aa_start = self._alt_data.variant_start_aa
                    while self._ref_seq[start] == self._alt_seq[start]:
                        start += 1
                        aa_start += 1
                    insertion = list(self._alt_seq[start:])
                    deletion = list(self._ref_seq[start:])
                    variants.append({"start": aa_start, "ins": insertion, "del": deletion})

                else:    # non-frameshifting delins or dup
                    # get size diff from diff in ref/alt lengths
                    start = self._alt_data.variant_start_aa - 1
                    aa_start = self._alt_data.variant_start_aa
                    delta = len(self._alt_seq) - len(self._ref_seq)
                    while self._ref_seq[start] == self._alt_seq[start]:
                        start += 1
                        aa_start += 1
                    offset = start + abs(delta)

                    if delta > 0:    # net insertion
                        insertion = list(self._alt_seq[start:offset])
                        deletion = []
                        ref_sub = self._ref_seq[start:]
                        alt_sub = self._alt_seq[offset:]
                    elif delta < 0:    # net deletion
                        insertion = []
                        deletion = list(self._ref_seq[start:offset])
                        ref_sub = self._ref_seq[offset:]
                        alt_sub = self._alt_seq[start:]
                    else:
                        insertion = []
                        deletion = []
                        ref_sub = self._ref_seq[start:]
                        alt_sub = self._alt_seq[start:]

                    # from start, get del/ins out to last difference
                    diff_indices = [i for i in xrange(len(ref_sub)) if ref_sub[i] != alt_sub[i]]
                    if diff_indices:
                        max_diff = diff_indices[-1] + 1
                        insertion.extend(list(alt_sub[:max_diff]))
                        deletion.extend(list(ref_sub[:max_diff]))

                    variants.append({"start": aa_start, "ins": insertion, "del": deletion})

            if DBG:
                print(variants)

        if self._is_ambiguous:
            var_ps = [self._create_variant('', '', '', '',
                                           acc=self._protein_accession,
                                           is_ambiguous=self._is_ambiguous)]
        elif len(self._alt_seq) == 0:
            var_ps = [self._create_variant('', '', '', '',
                                           acc=self._protein_accession,
                                           is_ambiguous=self._is_ambiguous,
                                           is_no_protein=True)]
        elif variants:
            var_ps = [self._convert_to_sequence_variants(x, self._protein_accession) for x in variants]
        else:    # ref = alt - "silent" hgvs change
            var_ps = [self._create_variant('', '', '', '', acc=self._protein_accession)]

        # TODO - handle multiple variants

        if len(var_ps) > 1:
            raise hgvs.exceptions.HGVSError("Got multiple AA variants - not supported")
        return var_ps[0]

    #
    # internal methods
    #

    def _convert_to_sequence_variants(self, variant, acc):
        """Convert AA variant to an hgvs representation

        :param variant: contains start, del, and ins
        :type variant: dict
        :param acc: protein accession
        :type acc: str
        :return hgvs string
        :rtype str
        """
        start = variant['start']
        insertion = ''.join(variant['ins'])
        deletion = ''.join(variant['del'])

        # defaults
        is_dup = False    # assume not dup
        fsext_len = None    # fs or ext length
        is_sub = False
        is_ext = False

        if start == 1:    # initial methionine is modified
            aa_start = aa_end = hgvs.location.AAPosition(base=start, aa=deletion)
            ref = ''
            alt = ''
            self._is_ambiguous = True    # side-effect

        if insertion and insertion.find("*") == 0:    # stop codon at variant position
            aa_start = aa_end = hgvs.location.AAPosition(base=start, aa=deletion[0])
            ref = ''
            alt = '*'
            is_sub = True

        elif start == len(self._ref_seq):    # extension
            if self._alt_seq[-1] == '*':
                fsext_len = len(insertion) - len(deletion)    # don't include the former stop codon
            else:
                fsext_len = '?'
            subst_at_stop_codon = insertion[0]

            aa_start = aa_end = hgvs.location.AAPosition(base=start, aa='*')
            ref = ''
            alt = subst_at_stop_codon
            is_ext = True

        elif self._is_frameshift:    # frameshift
            aa_start = aa_end = hgvs.location.AAPosition(base=start, aa=deletion[0])
            ref = ''

            try:
                fsext_len = str(insertion.index("*") + 1)    # start w/ 1st change; ends w/ * (inclusive)
            except ValueError:
                fsext_len = "?"

            alt = insertion[0]

        else:    # no frameshift - sub/delins/dup
            if len(insertion) == len(deletion) == 1:    # substitution
                aa_start = aa_end = hgvs.location.AAPosition(base=start, aa=deletion)
                ref = ''
                alt = insertion
                is_sub = True

            elif len(deletion) > 0:    # delins OR deletion OR stop codon at variant position
                ref = deletion
                end = start + len(deletion) - 1
                if len(insertion) > 0:    # delins
                    aa_start = hgvs.location.AAPosition(base=start, aa=deletion[0])
                    if end > start:
                        aa_end = hgvs.location.AAPosition(base=end, aa=deletion[-1])
                    else:
                        aa_end = aa_start
                    alt = insertion

                else:    # deletion OR stop codon at variant position
                    if len(deletion) + start == len(self._ref_seq):    # stop codon at variant position
                        aa_start = hgvs.location.AAPosition(base=start, aa=deletion[0])
                        aa_end = hgvs.location.AAPosition(base=start, aa=deletion[0])
                        ref = ''
                        alt = '*'
                        is_sub = True
                    else:    # deletion
                        aa_start = hgvs.location.AAPosition(base=start, aa=deletion[0])
                        if end > start:
                            aa_end = hgvs.location.AAPosition(base=end, aa=deletion[-1])
                        else:
                            aa_end = aa_start
                        alt = None

            elif len(deletion) == 0:    # insertion OR duplication OR extension

                is_dup, dup_start = self._check_if_ins_is_dup(start, insertion)

                if is_dup:    # duplication
                    dup_end = dup_start + len(insertion) - 1
                    aa_start = hgvs.location.AAPosition(base=dup_start, aa=insertion[0])
                    aa_end = hgvs.location.AAPosition(base=dup_end, aa=insertion[-1])
                    ref = alt = None

                else:    # insertion
                    start -= 1
                    end = start + 1

                    aa_start = hgvs.location.AAPosition(base=start, aa=self._ref_seq[start - 1])
                    aa_end = hgvs.location.AAPosition(base=end, aa=self._ref_seq[end - 1])
                    ref = None
                    alt = insertion

            else:    # should never get here
                raise ValueError("unexpected variant: {}".format(variant))

        var_p = self._create_variant(aa_start, aa_end, ref, alt,
                                     fsext_len=fsext_len,
                                     is_dup=is_dup,
                                     acc=acc,
                                     is_ambiguous=self._is_ambiguous,
                                     is_sub=is_sub,
                                     is_ext=is_ext)

        return var_p

    def _check_if_ins_is_dup(self, start, insertion):
        """Helper to identify an insertion as a duplicate

        :param start: 1-based insertion start
        :type start: int
        :param insertion: sequence
        :type insertion: str
        :return (is duplicate, variant start)
        :rtype (bool, int)
        """
        is_dup = False    # assume no
        variant_start = None

        dup_candidate_start = start - len(insertion) - 1
        dup_candidate = self._ref_seq[dup_candidate_start:dup_candidate_start + len(insertion)]
        if insertion == dup_candidate:
            is_dup = True
            variant_start = dup_candidate_start + 1

        return is_dup, variant_start

    def _create_variant(self, start, end, ref, alt,
                        fsext_len=None,
                        is_dup=False,
                        acc=None,
                        is_ambiguous=False,
                        is_sub=False,
                        is_ext=False,
                        is_no_protein=False):
        """Creates a SequenceVariant object"""
        interval = hgvs.location.Interval(start=start, end=end)
        # Note - order matters
        if is_no_protein:
            edit = '0'
        elif is_ambiguous:
            edit = '?'
        elif is_sub:
            edit = hgvs.edit.AASub(ref=ref, alt=alt)
        elif is_ext:
            edit = hgvs.edit.AAExt(ref=ref, alt=alt, aaterm='*', length=fsext_len)
        elif self._is_frameshift:
            edit = hgvs.edit.AAFs(ref=ref, alt=alt, length=fsext_len)
        elif is_dup:
            edit = hgvs.edit.Dup()
        elif ref == alt == '':
            edit = '='
        else:
            edit = hgvs.edit.AARefAlt(ref=ref, alt=alt)
        posedit = hgvs.posedit.PosEdit(interval, edit)
        if not (is_ambiguous and start == ''):
            posedit.uncertain = hgvs.global_config.mapping.inferred_p_is_uncertain
        var_p = hgvs.variant.SequenceVariant(acc, 'p', posedit)

        return var_p


## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
## 
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
## 
##     http://www.apache.org/licenses/LICENSE-2.0
## 
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
