"""simplified imports for the hgvs package

`hgvs.easy` simplifies using the hgvs package by providing a single
import path and objects that are instantiated with common defaults.

With default logging levels, using hgvs.easy is straightforward and
requires no code changes::

    >> from hgvs.easy import parser, projector
    >> var_g = parser.parse("NC_000017.11:g.43091687delC")
    >> projector.relevant_transcripts(var_g)
    ['NM_007294.3', 'NM_007297.3', 'NR_027676.1', 'NM_007298.3', 'NM_007299.3', 'NM_007300.3']
    >> projector.g_to_t(var_g, "NM_007294.3")
    SequenceVariant(ac=NM_007294.3, type=c, posedit=3844del)

`hgvs.easy` also introduces new functional forms for common methods.
For example::

    >> from hgvs.easy import parse, get_relevant_transcripts, g_to_t
    >> var_g = parse("NC_000017.11:g.43091687delC")
    >> get_relevant_transcripts(var_g)
    ['NM_007294.3', 'NM_007297.3', 'NR_027676.1', 'NM_007298.3', 'NM_007299.3', 'NM_007300.3']
    >> g_to_t(var_g, "NM_007294.3")
    SequenceVariant(ac=NM_007294.3, type=c, posedit=3844del)


NOTE: A consequence of making imports easy is a loss of
configurability by the caller.  The database connection is made with
no arguments (i.e., `connect()`), so it honors the UTA_DB_URL and
HGVS_SEQREPO_DIR environment variables but is otherwise not
configurable by the caller.

"""

from hgvs import __version__, global_config    # noqa: F401
from hgvs.assemblymapper import AssemblyMapper
from hgvs.dataproviders.uta import connect
from hgvs.normalizer import Normalizer
from hgvs.parser import Parser
from hgvs.validator import Validator
from hgvs.variantmapper import VariantMapper

# provide standard abbreviated, short, and long names for instances
hp = parser = hgvs_parser = Parser()
hdp = hgvs_data_provider = connect()
vm = variant_mapper = hgvs_variant_mapper = VariantMapper(hgvs_data_provider)
am37 = hgvs_assembly_mapper_37 = AssemblyMapper(hgvs_data_provider, assembly_name='GRCh37')
am38 = projector = hgvs_assembly_mapper_38 = AssemblyMapper(
    hgvs_data_provider, assembly_name='GRCh38')
hn = normalizer = hgvs_normalizer = Normalizer(hgvs_data_provider)
hv = validator = hgvs_validator = Validator(hgvs_data_provider)

# functionalized forms of common methods
parse = parser.parse
normalize = normalizer.normalize
validate = validator.validate
c_to_g = projector.c_to_g
c_to_n = projector.c_to_n
c_to_p = projector.c_to_p
g_to_c = projector.g_to_c
g_to_n = projector.g_to_n
g_to_t = projector.g_to_t
n_to_c = projector.n_to_c
n_to_g = projector.n_to_g
t_to_g = projector.t_to_g
t_to_p = projector.t_to_p
get_relevant_transcripts = am38.relevant_transcripts

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
