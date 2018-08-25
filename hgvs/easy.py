"""simplified imports for the hgvs package

`hgvs.easy` simplifies using the hgvs package by providing a single
import path and objects that are instantiated with common defaults.
Furthermore, objects are "lazy loaded" on first use: for example, the
grammar is parsed upon the first use of the parser, and a database
connection is made on the first use of any class that requires it.

With default logging levels, using hgvs.easy is straightforward and
requires no code changes::

    >> from hgvs.easy import parser, projector
    >> var_g = parser.parse("NC_000017.11:g.43091687delC")
    >> projector.relevant_transcripts(var_g)
    ['NM_007294.3', 'NM_007297.3', 'NR_027676.1', 'NM_007298.3', 'NM_007299.3', 'NM_007300.3']


Turning on debugging reveals what is happening behind the scenes::

    >> import logging
    >> logging.basicConfig(level="DEBUG")
    
    >> from hgvs.easy import parser, projector
    
    # Using the parser loads the grammar on demand
    >> var_g = parser.parse("NC_000017.11:g.43091687delC")
    DEBUG:hgvs.utils.lazywrapper:Lazy loading `parser = LazyWrapper(lambda: Parser())`...
    
    # Nested laziness: AssemblyMapper is wrapped, and its initialization
    # requires a hgvs dataprovider, which is also wrapped.  Using the
    # projector therefore connects to the db on demand, then creates the projector.
    >> projector.relevant_transcripts(var_g)
    DEBUG:hgvs.utils.lazywrapper:Lazy loading `projector = LazyWrapper(lambda: AssemblyMapper(hdp,...
    DEBUG:hgvs.utils.lazywrapper:Lazy loading `hdp = LazyWrapper(lambda: connect())` ...
    ['NM_007294.3', 'NM_007297.3', 'NR_027676.1', 'NM_007298.3', 'NM_007299.3', 'NM_007300.3']


NOTE: A consequence of making imports easy is a loss of
configurability by the caller.  The database connection is made with
no arguments (i.e., `connect()`), so it honors the UTA_DB_URL and
HGVS_SEQREPO_DIR environment variables but is otherwise not
configurable by the caller.

"""

from hgvs import __version__, global_config    # flake8: noqa
from hgvs.assemblymapper import AssemblyMapper
from hgvs.dataproviders.uta import connect
from hgvs.normalizer import Normalizer
from hgvs.parser import Parser
from hgvs.utils.lazywrapper import LazyWrapper
from hgvs.validator import Validator
from hgvs.variantmapper import VariantMapper


hp   = parser                  = LazyWrapper(Parser)
hdp  = hgvs_data_provider      = LazyWrapper(connect)
vm   = hgvs_variant_mapper     = VariantMapper(hgvs_data_provider)
am37 = hgvs_assembly_mapper_37 = LazyWrapper(lambda: AssemblyMapper(hgvs_data_provider, assembly_name='GRCh37'))
am38 = hgvs_assembly_mapper_38 = LazyWrapper(lambda: AssemblyMapper(hgvs_data_provider, assembly_name='GRCh38'))
hn   = hgvs_normalizer         = Normalizer(hgvs_data_provider)
hv   = hgvs_validator          = Validator(hgvs_data_provider)

projector = am38


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
