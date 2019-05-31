from __future__ import absolute_import, division, print_function, unicode_literals

import itertools
import re
import unittest

from hgvs.exceptions import HGVSDataNotAvailableError
import hgvs.dataproviders.uta
import hgvs.edit
import hgvs.parser
import hgvs.sequencevariant
import hgvs.variantmapper
from six.moves import map


class LineIterator(object):
    """iterate over a stream, keeping last line and # lines read"""

    def __init__(self, fh, skip_comments=False):
        self.fh = fh
        self.lines_read = None
        self.last_line = None
        self.skip_comments = skip_comments

    def __iter__(self):
        self.lines_read = 0
        for line in self.fh:
            self.lines_read += 1
            if self.skip_comments and line.startswith("#"):
                continue
            self.last_line = line
            yield line


class CrossChecker(object):
    """base class for testing a group of related variants.  Checks all reasonable combinations of
    g <-> t --> p

    """

    def crosscheck_variant_group(self, variants):
        """crosscheck a group of variants; returns None if successful, otherwise a message"""

        assert all(isinstance(v, hgvs.sequencevariant.SequenceVariant) for v in variants)

        variants = sorted(variants, key=lambda v: v.type)
        binned_variants = {k: [] for k in "cgmnrp"}
        binned_variants.update(
            {g: list(gi)
             for g, gi in itertools.groupby(variants, lambda v: v.type)})
        binned_variants["t"] = binned_variants["c"] + binned_variants["n"]

        assert len(binned_variants["g"]) == len(
            [v.ac for v in binned_variants["g"]]), "variants have multiple alignments"

        # g -> t: for each g., map to each transcript accession.
        for g_var in binned_variants["g"]:
            for t_var in binned_variants["t"]:
                try:
                    r = self.vm.g_to_t(g_var, t_var.ac)
                except HGVSDataNotAvailableError:
                    continue
                if t_var != r:
                    return "g_to_t({g_var},{t_var.ac}): got {r}; expected {t_var}".format(
                        g_var=g_var, t_var=t_var, r=r)

        # t -> g: for each t., map to each genomic accession
        for t_var in binned_variants["t"]:
            for g_var in binned_variants["g"]:
                try:
                    r = self.vm.t_to_g(t_var, g_var.ac)
                except HGVSDataNotAvailableError:
                    continue
                if g_var != r:
                    return "t_to_g({t_var},{g_var.ac}): got {r}; expected {g_var}".format(
                        g_var=g_var, t_var=t_var, r=r)

        # c -> p: for each c., map to a protein variant and check whether it's in result set

        if binned_variants["p"]:
            for c_var in binned_variants["c"]:
                try:
                    r = self.vm.c_to_p(c_var)
                except HGVSDataNotAvailableError:
                    continue
                r.posedit.uncertain = False
                if isinstance(r.posedit.edit, hgvs.edit.AAFs):
                    r.posedit.edit.length = None    # Clinvar doesn't have distance to Ter, so remove it
                if r not in binned_variants["p"]:
                    return "c_to_p({c_var}): got {r}; expected on of {p_vars}".format(
                        c_var=c_var, p_vars=" ".join(map(str, binned_variants["p"])), r=r)

        return None


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
