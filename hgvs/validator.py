# -*- coding: utf-8 -*-
"""implements validation of hgvs variants

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import hgvs
import hgvs.parser
import hgvs.edit
import hgvs.variantmapper
from hgvs.exceptions import HGVSInvalidVariantError, HGVSUnsupportedOperationError
from hgvs.enums import ValidationLevel, Datum

SEQ_ERROR_MSG = "Variant reference ({var_ref_seq}) does not agree with reference sequence ({ref_seq})"
CDS_BOUND_ERROR_MSG = "Variant coordinate is out of the bound of CDS region (CDS length : {cds_length})"

BASE_OFFSET_COORD_TYPES = "cnr"
SIMPLE_COORD_TYPES = "gmp"


class Validator(object):
    """invoke intrinsic and extrinsic validation"""

    def __init__(self, hdp, strict=hgvs.global_config.validator.strict):
        self.strict = strict
        self._ivr = IntrinsicValidator(strict)
        self._evr = ExtrinsicValidator(hdp, strict)

    def validate(self, var, strict=None):
        if strict is None:
            strict = self.strict
        return self._ivr.validate(var, strict) and self._evr.validate(var, strict)


class IntrinsicValidator(object):
    """
    Attempts to determine if the HGVS name is internally consistent

    """

    def __init__(self, strict=hgvs.global_config.validator.strict):
        self.strict = strict

    def validate(self, var, strict=None):
        assert isinstance(var, hgvs.sequencevariant.SequenceVariant
                          ), "variant must be a parsed HGVS sequence variant object"
        if strict is None:
            strict = self.strict
        fail_level = ValidationLevel.WARNING if strict else ValidationLevel.ERROR
        (res, msg) = var.validate()
        if res >= fail_level:
            raise HGVSInvalidVariantError(msg)
        return True


class ExtrinsicValidator():
    """
    Attempts to determine if the HGVS name validates against external data sources
    """

    def __init__(self, hdp, strict=hgvs.global_config.validator.strict):
        self.strict = strict
        self.hdp = hdp
        self.vm = hgvs.variantmapper.VariantMapper(self.hdp, prevalidation_level=None)

    def validate(self, var, strict=None):
        assert isinstance(var, hgvs.sequencevariant.SequenceVariant
                          ), "variant must be a parsed HGVS sequence variant object"
        if strict is None: strict = self.strict
        fail_level = ValidationLevel.WARNING if strict else ValidationLevel.ERROR
        (res, msg) = self._ref_is_valid(var)
        if res >= fail_level:
            raise HGVSInvalidVariantError(msg)
        else:
            (res, msg) = self._c_within_cds_bound(var)
            if res >= fail_level:
                raise HGVSInvalidVariantError(msg)
        return True

    def _ref_is_valid(self, var):
        # use reference sequence of original variant, even if later converted (eg c_to_n)
        if (var.type in BASE_OFFSET_COORD_TYPES and var.posedit.pos is not None
                and (var.posedit.pos.start.offset != 0 or var.posedit.pos.end.offset != 0)):
            return (ValidationLevel.WARNING,
                    "Cannot validate sequence of an intronic variant ({})".format(str(var)))

        ref_checks = []
        if var.type == 'p':
            if not var.posedit or not var.posedit.pos or not var.posedit.pos.start or not var.posedit.pos.end:
                return (ValidationLevel.VALID, None)
            ref_checks.append((var.ac, var.posedit.pos.start.pos, var.posedit.pos.start.pos,
                               var.posedit.pos.start.aa))
            if var.posedit.pos.start.pos != var.posedit.pos.end.pos:
                ref_checks.append((var.ac, var.posedit.pos.end.pos, var.posedit.pos.end.pos,
                                   var.posedit.pos.end.aa))
        else:
            var_ref_seq = getattr(var.posedit.edit, "ref", None) or None
            var_x = self.vm.c_to_n(var) if var.type == "c" else var
            ref_checks.append((var_x.ac, var_x.posedit.pos.start.base, var_x.posedit.pos.end.base,
                               var_ref_seq))

        for ac, var_ref_start, var_ref_end, var_ref_seq in ref_checks:
            if var_ref_start is None or var_ref_end is None or not var_ref_seq:
                continue

            # ref_seq is digit, as in "del6"
            try:
                int(var_ref_seq)
                continue
            except ValueError:
                pass

            ref_seq = self.hdp.get_seq(ac, var_ref_start - 1, var_ref_end)
            if ref_seq != var_ref_seq:
                return (ValidationLevel.ERROR, str(var) + ": " + SEQ_ERROR_MSG.format(
                    ref_seq=ref_seq, var_ref_seq=var_ref_seq))

        return (ValidationLevel.VALID, None)

    def _c_within_cds_bound(self, var):
        if var.type != 'c':
            return (ValidationLevel.VALID, None)
        tx_info = self.hdp.get_tx_identity_info(var.ac)
        if tx_info is None:
            return (ValidationLevel.WARNING,
                    "No transcript data for accession: {ac}".format(ac=var.ac))
        cds_length = tx_info["cds_end_i"] - tx_info["cds_start_i"]
        if var.posedit.pos.start.datum == Datum.CDS_START and var.posedit.pos.start.base > cds_length:
            return (ValidationLevel.ERROR, CDS_BOUND_ERROR_MSG.format(cds_length=cds_length))
        if var.posedit.pos.end.datum == Datum.CDS_START and var.posedit.pos.end.base > cds_length:
            return (ValidationLevel.ERROR, CDS_BOUND_ERROR_MSG.format(cds_length=cds_length))
        return (ValidationLevel.VALID, None)


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
