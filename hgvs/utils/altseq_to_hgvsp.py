#
# Utility class for creating an hgvsp SequenceVariant object,
# given a transcript with variants applied.
# Used in hgvsc to hgvsp conversion.
#
# The "frameshift_start" param is used to identify to the code
# where any "terminal" frameshift should occur.
# Consider the case where there is a single base deletion, but the resulting
# frameshift leaves an island of amino acids that coincidentally align to the reference.
# This can lead to splitting a single frameshift call into [indel+match+downstream frameshift]
# This param forces the comparison code to mark the rest of the sequence as part of a frameshift
# and terminates the comparison.
#
import collections
import difflib

import hgvs.edit
import hgvs.exceptions
import hgvs.location
import hgvs.posedit
import hgvs.utils
import hgvs.variant

DBG = False

class AltSeqToHgvsp(object):

    def __init__(self, ref_seq, alt_seq, protein_accession, frameshift_start=None, is_substitution=False,
                 is_ambiguous=True):
        """Constructor

        :param ref_seq: protein reference sequence
        :type str
        :param alt_seq: protein alternate sequence
        :type str
        :param frameshift_start: index n protein (1-based) where a frameshift should affect all downstream bases
        :type int
        :param is_substitution: do a simple compare if you expect a single AA difference
        :type bool
        :param is_ambiguous: cannot predict protein
        type bool
        """
        self._ref_seq = ref_seq
        self._alt_seq = alt_seq
        self._frameshift_start = frameshift_start
        self._protein_accession = protein_accession
        self._is_substitution = is_substitution
        self._is_ambiguous = is_ambiguous

        if DBG:
            print len(ref_seq), len(alt_seq), frameshift_start, protein_accession
            print ref_seq
            print alt_seq

    def build_hgvsp(self):
        """Compare two amino acid sequences; generate an hgvs tag from the output

        :return list of variants in sequence order
        :type list of dict
        """

        indel_map = {'+':'ins', '-':'del'}
        variants = []

        if not self._is_ambiguous:
            # for a simple substitution, just do a straight AA-by-AA comparison
            do_difflib = True   # assume complex case
            if self._is_substitution and len(self._ref_seq) == len(self._alt_seq):
                diff_pos = [(i, self._ref_seq[i], self._alt_seq[i]) for i in xrange(len(self._ref_seq))
                            if  self._ref_seq[i] != self._alt_seq[i]]
                if len(diff_pos) == 1:
                    (start, deletion, insertion) = diff_pos[0]
                    do_difflib = False
                    variants.append({"start": start + 1, "ins": insertion, "del": deletion})

            # simple comparison didn't work or isn't appropriate - resort to difflib
            if do_difflib:
                diff_list = difflib.Differ().compare(self._ref_seq, self._alt_seq)

                # walk through the difflib string list and summarize variants as they are encountered
                ref_index = 1
                alt_index = 1
                in_variant = False
                if DBG:
                    item_dbg = []
                for item in diff_list:
                    if DBG:
                        item_dbg.append(item)
                    change = item[0]
                    if change == ' ':       # a match
                        if in_variant:      # process the variant just exited
                            in_variant = False
                            variants.append(current_var)
                        ref_index += 1
                        alt_index += 1
                    else:                   # a mismatch
                        base = item[-1]
                        if not in_variant:  # start a new variant
                            in_variant = True
                            current_var = collections.defaultdict(list)
                            current_var['start'] = ref_index

                            if self._frameshift_start and current_var['start'] >= self._frameshift_start:
                                current_var = self._force_variant_to_frameshift(current_var, ref_index, alt_index)
                                variants.append(current_var)
                                in_variant = False
                                break

                        current_var[indel_map[change]].append(base)
                        if change == '-':
                            ref_index += 1
                        elif change == '+':
                            alt_index += 1


                if in_variant:  # implies end of a frameshift; need to add this last variant to the list
                    variants.append(current_var)

                if DBG:
                    print item_dbg
                    print variants

        if self._is_ambiguous:
            var_ps = [self._create_variant('', '', '', '', acc=self._protein_accession,
                                           is_ambiguous=self._is_ambiguous)]
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

    def _force_variant_to_frameshift(self, current_variant, ref_index, alt_index):
        """Make this variant a terminal frameshift

        :param current_variant: representation of variant to describe as a frameshift
        :type dict
        :param ref_index: reference sequence index (1-based) where frameshift starts
        :type int
        :param alt_index:  alt sequence index (1-based) where frameshift starts
        :type int
        :return updated current_variant
        :type dict
        """
        current_variant['del'] = list(self._ref_seq[ref_index - 1:])
        current_variant['ins'] = list(self._alt_seq[alt_index - 1:])
        return current_variant

    def _convert_to_sequence_variants(self, variant, acc):
        """Convert AA variant to an hgvs representation

        :param variant: contains start, del, and ins
        :type dict
        :param acc: protein accession
        :type str
        :return hgvs string
        """
        start = variant['start']
        insertion = ''.join(variant['ins'])
        deletion = ''.join(variant['del'])

        is_frameshift = insertion and deletion and deletion[-1] == '*'

        # defaults
        is_dup = False  # assume not dup
        fs = None

        if is_frameshift:                                               # frameshift
            aa_start = aa_end = hgvs.location.AAPosition(base=start, aa=deletion[0])
            ref = ''

            try:
                new_stop = str(insertion.index("*") + 1)    # start w/ 1st change; ends w/ * (inclusive)
            except ValueError:
                new_stop = "?"

            if new_stop != "1":
                alt = insertion[0]
                fs = 'fs*{}'.format(new_stop)

            else:   # frameshift introduced stop codon at variant position
                alt = '*'
        elif start == 1:                                          # initial methionine is modified
                aa_start = aa_end = hgvs.location.AAPosition(base=start, aa=deletion)
                ref = ''
                alt = ''
                self._is_ambiguous = True   # side-effect

        else:                                                           # no frameshift
            if len(insertion) == len(deletion) == 1:                    # substitution
                aa_start = aa_end = hgvs.location.AAPosition(base=start, aa=deletion)
                ref = ''
                alt = insertion

            elif len(deletion) > 0:                                     # delins OR deletion
                ref = deletion


                end = start + len(deletion) - 1
                if len(insertion) > 0:                                  # delins
                    aa_start = hgvs.location.AAPosition(base=start, aa=deletion[0])
                    if end > start:
                        aa_end =  hgvs.location.AAPosition(base=end, aa=deletion[-1])
                    else:
                        aa_end = aa_start
                    alt = insertion

                else:                                                   # deletion OR stop codon at variant position
                    if len(deletion) + start == len(self._ref_seq):     # stop codon at variant position
                        aa_start = hgvs.location.AAPosition(base=start, aa=deletion[0])
                        aa_end = hgvs.location.AAPosition(base=start, aa=deletion[0])
                        ref = ''
                        alt = '*'
                    else:                                               # deletion
                        # get c-terminal representation
                        c_term_start = self._get_shift_match(self._ref_seq, start - 1, deletion) + 1
                        c_term_end = c_term_start + (end - start)

                        aa_start = hgvs.location.AAPosition(base=c_term_start, aa=deletion[0])
                        if end > start:
                            aa_end =  hgvs.location.AAPosition(base=c_term_end, aa=deletion[-1])
                        else:
                            aa_end = aa_start
                        alt = None

            elif len(deletion) == 0:                                    # insertion OR duplication OR extension

                is_dup, dup_start = self._check_if_ins_is_dup(start, insertion)

                if is_dup:                                              # duplication
                    # get c-terminal
                    dup_start = self._get_shift_match(self._ref_seq, dup_start -1, insertion) + 1

                    dup_end = dup_start + len(insertion) - 1
                    aa_start = hgvs.location.AAPosition(base=dup_start, aa=insertion[0])
                    aa_end =  hgvs.location.AAPosition(base=dup_end, aa=insertion[-1])
                    ref = alt = None

                else:                                                   # insertion OR extension
                    if start == len(self._ref_seq):                     # extension
                        len_ext = len(insertion) # don't include the former stop codon
                        subst_at_stop_codon = insertion[-1]

                        aa_start = aa_end = hgvs.location.AAPosition(base=start, aa='*')
                        ref = ''
                        alt = subst_at_stop_codon
                        fs ='ext*{}'.format(len_ext)

                    else:                                               # insertion
                        start = start - 1
                        end = start + 1

                        aa_start = hgvs.location.AAPosition(base=start, aa=self._ref_seq[start - 1])
                        aa_end =  hgvs.location.AAPosition(base=end, aa=self._ref_seq[end - 1])
                        ref = None
                        alt = insertion

            else: # should never get here
                raise ValueError("unexpected variant: {}".format(variant))

        var_p = self._create_variant(aa_start, aa_end, ref, alt, fs, is_dup, acc, is_ambiguous=self._is_ambiguous)

        return var_p


    def _check_if_ins_is_dup(self, start, insertion):
        """Helper to identify an insertion as a duplicate"""
        is_dup = False  # assume no
        variant_start = None

        if len(insertion) == 1:  # since difflib may match 1-AA dups before or after, need to check both
            ref_prev = self._ref_seq[start - 2]
            ref_next = self._ref_seq[start - 1]
            if insertion == ref_prev or insertion == ref_next:
                is_dup = True
                variant_start = start - 1 if insertion == ref_prev else start
        else:
            # figure out the candidate ref positions that should match if the insertion is a dup
            dup_candidate_start = start - len(insertion) - 1
            dup_candidate_end = dup_candidate_start + len(insertion)
            dup_candidate = self._ref_seq[dup_candidate_start: dup_candidate_end]

            if insertion == dup_candidate:
                is_dup = True
                variant_start = dup_candidate_start + 1

        return is_dup, variant_start

    def _get_shift_match(self, seq, start, match_seq, is_rev=False):
        """Helper to identify a terminal match of a string in a sequence"""
        result = start
        sign = -1 if is_rev else 1
        result = self._find_shift(seq, sign*len(match_seq), match_seq, result)
        result = self._find_shift(seq, sign*1, match_seq, result)
        return result

    def _find_shift(self, seq, increment, match_seq, init_result):
        """Helper for _get_shift_match - scans over an interval for the most extreme match"""
        match_failed = False
        cur_start = init_result
        while not match_failed:
            result = cur_start
            cur_start += increment
            match_failed = seq[cur_start:cur_start + len(match_seq)] != match_seq
        return result

    def _create_variant(self, start, end, ref, alt, fs=None, is_dup=False, acc=None, is_ambiguous=False):
        """Creates a SequenceVariant object"""
        interval = hgvs.location.Interval(start=start, end=end)
        if is_ambiguous:
            edit = hgvs.edit.AASpecial(status='?')
        elif is_dup:
            edit = hgvs.edit.Dup()
        elif ref == alt == '':
            edit = hgvs.edit.AASpecial(status='=')
        else:
            edit = hgvs.edit.AARefAlt(ref=ref, alt=alt, fs=fs)
        posedit = hgvs.posedit.PosEdit(interval, edit)
        var_p = hgvs.variant.SequenceVariant(acc, 'p', posedit)

        return var_p

