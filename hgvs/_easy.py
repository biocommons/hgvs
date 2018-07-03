"""hgvs._easy -- one-shot imports for hgvs

WARNING: This is for internal use only. It is subject to removal without warning.

Example:
>>> import hgvs._easy as hgvs
>>> hgvs.__version__
>>> hp = hgvs.Parser()

"""

from hgvs import __version__    # flake8: noqa
from hgvs.assemblymapper import AssemblyMapper
from hgvs.dataproviders.uta import connect
from hgvs.normalizer import Normalizer
from hgvs.parser import Parser
from hgvs.validator import Validator
from hgvs.variantmapper import VariantMapper


if False:
    hp = hgvs_parser = Parser()
    hdp = hgvs_data_provider = connect()
    vm = hgvs_variant_mapper = VariantMapper(hgvs_data_provider)
    am37 = hgvs_assembly_mapper_37 = AssemblyMapper(hgvs_data_provider, assembly_name='GRCh37')
    am38 = hgvs_assembly_mapper_38 = AssemblyMapper(hgvs_data_provider, assembly_name='GRCh38')
    hn = hgvs_normalizer = Normalizer(hgvs_data_provider)
    hv = hgvs_validator = Validator(hgvs_data_provider)


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
