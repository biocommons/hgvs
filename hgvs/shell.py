# -*- coding: utf-8 -*-
"""start IPython shell with hgvs initialized. Intended to be used for
experimenting, debugging, and generating bug reports."""

from __future__ import absolute_import, division, print_function, unicode_literals

import logging

import IPython

header_string = """############################################################################
hgvs-shell -- interactive hgvs
hgvs version: {v}
data provider url: {hdp.url}
schema_version: {sv}
data_version: {dv}
sequences source: {hdp.seqfetcher.source}

The following variables are defined:
* global_config
* hp, parser, hgvs_parser -- Parser instance
* hdp, hgvs_data_provider -- UTA data provider instance
* vm, variant_mapper, hgvs_variant_mapper -- VariantMapper instance
* am37, hgvs_assembly_mapper_37 -- GRCh37 Assembly Mapper instance
* am38, projector, hgvs_assembly_mapper_38 -- GRCh38 Assembly Mapper instances
* hn, normalizer, hgvs_normalizer -- Normalizer instance
* hv, validator, hgvs_validator) -- Validator instance

The following functions are available:
  * parse, normalize, validate
  * g_to_c, g_to_n, g_to_t,
  * c_to_g, c_to_n, c_to_p,
  * n_to_c, n_to_g,
  * t_to_g,
  * get_relevant_transcripts

When submitting bug reports, include the version header shown above
and use these variables/variable names whenever possible.

"""


def shell():
    logging.basicConfig(level=logging.WARNING)

    from hgvs.easy import (    # noqa: F401
        __version__, global_config,

    # instances
        hp, parser, hgvs_parser, hdp, hgvs_data_provider, vm, variant_mapper, hgvs_variant_mapper,
        am37, hgvs_assembly_mapper_37, am38, projector, hgvs_assembly_mapper_38, hn, normalizer,
        hgvs_normalizer, hv, validator, hgvs_validator,

    # functionalized methods
        parse, normalize, validate, g_to_c, g_to_n, g_to_t, c_to_g, c_to_n, c_to_p, n_to_c, n_to_g,
        t_to_g, t_to_p, get_relevant_transcripts)

    from hgvs.utils.context import variant_context_w_alignment    # noqa

    IPython.embed(
        header=header_string.format(
            v=__version__,
            hdp=hdp,
            sv=hdp.schema_version(),
            dv=hdp.data_version(),
        ))


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
