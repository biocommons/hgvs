#
# Compares two protein sequences
#

import collections
import difflib
import pprint

import hgvs.edit as edit
import hgvs.location as location
import hgvs.translators.utils.aminoacidutils as aminoacidutils


DBG = True

class ProteinComparer(object):

    def __init__(self, ref_seq, alt_seq, frameshift_start = -1, force_frameshift=True):
        """Constructor

        The normal difflib comparison may identify islands of matches in a frameshift
        which will end up being reported as multiple variants.

        The force_frameshift param is used to tell the app to treat that frameshift as the last variant,
        rather than just use the results of the diff.

        :param ref_seq: protein reference sequence
        :type str
        :param alt_seq: protein alternate sequence
        :type str
        :param frameshift_start: index n protein (1-based) where a frameshift should affect all downstream bases
        :type int
        :param force_frameshift: treat a variant which would shift all subsequent bases as the terminal frameshift
        :type bool
        """
        self._ref_seq = ref_seq
        self._alt_seq = alt_seq
        self._frameshift_start = frameshift_start
        self._force_frameshift = force_frameshift

    def compare(self):
        """Compare two amino acid sequences; generate an hgvs tag from the output

        :return list of variants in sequence order
        :type list of dict
        """

        indel_map = {'+':'ins', '-':'del'}

        diff_list = difflib.Differ().compare(self._ref_seq, self._alt_seq)

        # initial check for frameshift
        ref_stop = len(self._ref_seq) - 1
        alt_stop = self._alt_seq.find("*")
        has_frameshift = alt_stop != -1 and alt_stop != ref_stop


        # walk through the difflib string list and summarize variants as they are encountered
        ref_index = 1
        alt_index = 1
        in_variant = False
        variants = []
        if DBG:
            debug_item = []
        for item in diff_list:
            if DBG:
                debug_item.append(item)
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

                    if self._force_frameshift:
                        if self._frameshift_start != -1 and current_var['start'] >= self._frameshift_start:
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
            print(debug_item)
            pprint.pprint(variants)

        #sequence_variants = [self._convert_to_sequence_variants(x) for x in variants]

        return variants

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

    # def _convert_to_sequence_variants(self, variant):
    #     """Convert AA variant (as a dictionary) to SequenceVariant
    #     :param variant: contains start, del, and ins
    #     :type dict
    #     :return SequenceVariant
    #     """
    #     start = variant['start']
    #     insertion = variant['ins']
    #     deletion = variant['del']
    #
    #     is_frameshift = insertion and deletion and deletion[-1] == '*'
    #
    #     # frameshift
    #     if is_frameshift:
    #         delstart3 = aminoacidutils.convert_AA_1_to_3(deletion[0])
    #         insstart3 = aminoacidutils.convert_AA_1_to_3(insertion[0])
    #         try:
    #             new_stop = str(insertion.index("*") + 1)    # start w/ 1st change; ends w/ * (inclusive)
    #         except ValueError:
    #             new_stop = "?"
    #
    #         if new_stop != "1":
    #             tag = "{}{}{}fs*{}".format(delstart3, start, insstart3, new_stop)
    #         else:   # frameshift introduced stop codon at variant position
    #             tag = "{}{}*".format(delstart3, start)
    #
    #     # no frameshift
    #     else:
    #
    #         deletion3 = aminoacidutils.convert_AA_1_to_3(deletion)
    #         insertion3 = aminoacidutils.convert_AA_1_to_3(insertion)
    #
    #         if len(insertion) == len(deletion) == 1:        # substitution
    #             tag = '{}{}{}'.format(deletion3, start, insertion3)
    #             ref_alt = edit.RefAlt(ref=None, alt=insertion3)
    #             aapos_start = location.AAPosition(pos=)
    #
    #         elif len(deletion) > 0:
    #             end = start + len(deletion) - 1                         # delins OR deletion
    #             delstart3 = aminoacidutils.convert_AA_1_to_3(deletion[0])
    #             delend3 = aminoacidutils.convert_AA_1_to_3(deletion[-1])
    #
    #             if len(insertion) > 0:                                  # delins
    #                 tag = '{}{}_{}{}delins{}'.format(delstart3, start, delend3, end, insertion3) if end > start \
    #                     else '{}{}delins{}'.format(deletion3, start, insertion3)
    #             else:                                                   # deletion
    #                 if len(deletion) + start == len(self._ref_seq): # stop codon introduced at variant position
    #                     tag = '{}{}*'.format(delstart3, start)
    #                 else:
    #                     tag = '{}{}_{}{}del'.format(delstart3, start, delend3, end) if end > start \
    #                         else '{}{}del'.format(deletion3, start)
    #
    #         elif len(deletion) == 0:                                    # insertion OR duplication
    #
    #             is_dup, variant_start = self._check_if_ins_is_dup(start, insertion, self._ref_seq)
    #
    #             if is_dup:                                              # duplication
    #                 dupstart3 = aminoacidutils.convert_AA_1_to_3(insertion[0])
    #                 dupend3 = aminoacidutils.convert_AA_1_to_3(insertion[-1])
    #                 tag = '{}{}_{}{}dup'.format(dupstart3, variant_start, dupend3, variant_start + len(insertion) - 1) \
    #                     if len(insertion) > 1 else '{}{}dup'.format(dupstart3, variant_start)
    #
    #
    #             else:                                                   # insertion
    #                 if start == len(self._ref_seq):                           # extension
    #                     len_ext = len(insertion) - 1  # don't include the former stop codon
    #                     subst_at_stop_codon = aminoacidutils.convert_AA_1_to_3(insertion[-1])
    #                     tag = '*{}{}ext*{}'.format(start, subst_at_stop_codon, len_ext)
    #                 else:                                               # insertion
    #                     ref_start3 = aminoacidutils.convert_AA_1_to_3(self._ref_seq[start - 1])
    #                     ref_end3 = aminoacidutils.convert_AA_1_to_3(self._ref_seq[start])
    #                     tag = '{}{}_{}{}ins{}'.format(ref_start3, start, ref_end3, start + 1, insertion3)
    #
    #         else: # should never get here
    #             raise ValueError("unexpected variant: {}".format(variant))
    #
    #     return tag
    #
    #
    # def _check_if_ins_is_dup(self, start, insertion, ref_seq):
    #     """Helper to identify an insertion as a duplicate
    #     :param start 1-based start position in reference sequence
    #     :type int
    #     :param insertion amino acid insertion sequence
    #     :type list
    #     :param ref_seq amino acid reference sequence
    #     :type str
    #     :return (is duplicate [boolean], variant start [-1 if not a duplicate)
    #     :type tuple
    #     """
    #
    #     is_dup = False  # assume no
    #     variant_start = -1
    #
    #     if len(insertion) == 1:  # since difflib may match single AA dups before or after, need to check both
    #         ref_prev = [ref_seq[start - 2]]
    #         ref_next = [ref_seq[start - 1]]
    #         if DBG:
    #             print start
    #             print insertion
    #             print ref_prev, ref_next
    #         if insertion == ref_prev or insertion == ref_next:
    #             is_dup = True
    #             variant_start = start - 1 if insertion == ref_prev else start
    #
    #     else:
    #         # figure out the candidate ref positions that should match if the insertion is a dup
    #         dup_candidate_start = start - len(insertion) - 1
    #         dup_candidate_end = dup_candidate_start + len(insertion)
    #         dup_candidate = list(ref_seq[dup_candidate_start: dup_candidate_end])
    #
    #         if DBG:
    #             print dup_candidate_start, dup_candidate_end
    #             pprint.pprint(insertion)
    #             pprint.pprint(dup_candidate)
    #
    #         if insertion == dup_candidate:
    #             is_dup = True
    #             variant_start = dup_candidate_start + 1
    #
    #     return is_dup, variant_start



def main():
    pass


if __name__ == "__main__":
    main()

