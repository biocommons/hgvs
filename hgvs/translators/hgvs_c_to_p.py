#
# Translate hgvsc tag to hgvsp tag
#
import pprint

from Bio.Seq import Seq

import hgvs.edit
import hgvs.parser
import hgvs.translators.utils.aminoacidutils as aminoacidutils
import hgvs.translators.utils.proteincomparer as proteincomparer
import hgvs.translators.utils.variantinserter as variantinserter

DBG = True  # DEBUG flag

class HgvsCToP(object):


    def __init__(self, datasource):
        """Constructor

        :param datasource
        :type base_input_source
        """
        self.datasource = datasource
        self.parser = hgvs.parser.Parser()

    def convert(self, hgvsc):
        """Convert hgvsc tag to hgvsp tag

        :param hgvsc: hgvsc tag
        :type str
        """

        # translate hgvsc
        var = self.parser.hgvs_variant(hgvsc)
        if DBG:
            print "START DEBUG"
            print hgvsc

        # get transcript from datasource & convert to AA
        transcript_data = self.datasource.get_sequence(var.seqref)

        transcript_dna_cds = \
            Seq(transcript_data['transcript_sequence'][transcript_data['cds_start'] - 1:transcript_data['cds_stop']])
        transcript_data["AA_sequence"] = str(transcript_dna_cds.translate())

        # incorporate DNA change into sequence(s) and translate
        inserter = variantinserter.VariantInserter(var, transcript_data)
        variant_data = inserter.insert_variant()

        if DBG:
            pprint.pprint(transcript_data)
            pprint.pprint(variant_data[0])

        # perform comparison to get hgvs tag
        comparer = proteincomparer.ProteinComparer(transcript_data['AA_sequence'],
                                                variant_data[0]['AA_sequence'], variant_data[0]['frameshift_start'])
        vars = comparer.compare()
        hgvsp = self._convert_all_to_hgvsp(vars, transcript_data['AA_sequence'], variant_data[0]['AA_sequence'],
                                           variant_data[0]['protein_accession'])

        if DBG:
            print "END DEBUG"

        return hgvsp

    #
    # internal methods
    #

    def _convert_all_to_hgvsp(self, variants, ref_seq, alt_seq, accession):


        # convert variant to hgvs p. format
        # p.[(<variant1>, <variant2,...>)]
        hgvsps = [self._convert_AA_variant_to_hgvs(variant, ref_seq, alt_seq) for variant in variants]

        if len(hgvsps) > 1:
            result = "{}:p.[({})]".format(accession, ';'.join(hgvsps))
        elif len(hgvsps) == 1:
            result = "{}:p.({})".format(accession, hgvsps[0])
        else:
            result = "{}:p.(=)".format(accession)

        if DBG:
            print hgvsps
            print result

        return result


    def _convert_AA_variant_to_hgvs(self, variant, ref_seq, alt_seq):
        """Convert AA variant to an hgvs representation
        :param variant: contains start, del, and ins
        :type dict
        :param ref_seq: AA reference
        :type string
        :param alt_seq: AA reference with variants
        :type string
        :return hgvs string
        """
        start = variant['start']
        insertion = variant['ins']
        deletion = variant['del']

        is_frameshift = insertion and deletion and deletion[-1] == '*'

        # frameshift
        if is_frameshift:
            delstart3 = aminoacidutils.convert_AA_1_to_3(deletion[0])
            insstart3 = aminoacidutils.convert_AA_1_to_3(insertion[0])
            try:
                new_stop = str(insertion.index("*") + 1)    # start w/ 1st change; ends w/ * (inclusive)
            except ValueError:
                new_stop = "?"

            if new_stop != "1":
                tag = "{}{}{}fs*{}".format(delstart3, start, insstart3, new_stop)
            else:   # frameshift introduced stop codon at variant position
                tag = "{}{}*".format(delstart3, start)

        # no frameshift
        else:

            deletion3 = aminoacidutils.convert_AA_1_to_3(deletion)
            insertion3 = aminoacidutils.convert_AA_1_to_3(insertion)

            if len(insertion) == len(deletion) == 1:        # substitution
                tag = '{}{}{}'.format(deletion3, start, insertion3)

            elif len(deletion) > 0:
                end = start + len(deletion) - 1                         # delins OR deletion
                delstart3 = aminoacidutils.convert_AA_1_to_3(deletion[0])
                delend3 = aminoacidutils.convert_AA_1_to_3(deletion[-1])

                if len(insertion) > 0:                                  # delins
                    tag = '{}{}_{}{}delins{}'.format(delstart3, start, delend3, end, insertion3) if end > start \
                        else '{}{}delins{}'.format(deletion3, start, insertion3)
                else:                                                   # deletion
                    if len(deletion) + start == len(ref_seq): # stop codon introduced at variant position
                        tag = '{}{}*'.format(delstart3, start)
                    else:
                        tag = '{}{}_{}{}del'.format(delstart3, start, delend3, end) if end > start \
                            else '{}{}del'.format(deletion3, start)

            elif len(deletion) == 0:                                    # insertion OR duplication

                is_dup, variant_start = self._check_if_ins_is_dup(start, insertion, ref_seq)

                if is_dup:                                              # duplication
                    dupstart3 = aminoacidutils.convert_AA_1_to_3(insertion[0])
                    dupend3 = aminoacidutils.convert_AA_1_to_3(insertion[-1])
                    tag = '{}{}_{}{}dup'.format(dupstart3, variant_start, dupend3, variant_start + len(insertion) - 1) \
                        if len(insertion) > 1 else '{}{}dup'.format(dupstart3, variant_start)


                else:                                                   # insertion
                    if start == len(ref_seq):                           # extension
                        len_ext = len(insertion) - 1  # don't include the former stop codon
                        subst_at_stop_codon = aminoacidutils.convert_AA_1_to_3(insertion[-1])
                        tag = '*{}{}ext*{}'.format(start, subst_at_stop_codon, len_ext)
                    else:                                               # insertion
                        ref_start3 = aminoacidutils.convert_AA_1_to_3(ref_seq[start - 1])
                        ref_end3 = aminoacidutils.convert_AA_1_to_3(ref_seq[start])
                        tag = '{}{}_{}{}ins{}'.format(ref_start3, start, ref_end3, start + 1, insertion3)

            else: # should never get here
                raise ValueError("unexpected variant: {}".format(variant))

        return tag


    def _check_if_ins_is_dup(self, start, insertion, ref_seq):
        """Helper to identify an insertion as a duplicate
        """

        is_dup = False  # assume no
        variant_start = -1

        if len(insertion) == 1:  # since difflib may match 1-AA dups before or after, need to check both
            ref_prev = [ref_seq[start - 2]]
            ref_next = [ref_seq[start - 1]]
            if DBG:
                print start
                print insertion
                print ref_prev, ref_next
            if insertion == ref_prev or insertion == ref_next:
                is_dup = True
                variant_start = start - 1 if insertion == ref_prev else start

        else:
            # figure out the candidate ref positions that should match if the insertion is a dup
            dup_candidate_start = start - len(insertion) - 1
            dup_candidate_end = dup_candidate_start + len(insertion)
            dup_candidate = list(ref_seq[dup_candidate_start: dup_candidate_end])

            if DBG:
                print dup_candidate_start, dup_candidate_end
                pprint.pprint(insertion)
                pprint.pprint(dup_candidate)

            if insertion == dup_candidate:
                is_dup = True
                variant_start = dup_candidate_start + 1

        return is_dup, variant_start



def main():
    pass


if __name__ == "__main__":
    main()

