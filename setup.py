from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

import re


with open('doc/description.txt') as f:
    long_description = f.read()


def version_handler(mgr, options):
    version = mgr.get_current_version()
    if version.endswith('dev'):
        version += '-' + mgr._invoke('log', '-l1', '-r.', '--template', '{node|short}').strip()
    elif re.match('^\d+\.\d+$', version):
        # StrictVersion considers x.y == x.y.0 and drops the .0 from a
        # repo tag.  Add it back and ensure that it's really a tag for
        # our parent.
        version += '.0'
        assert version in mgr.get_parent_tags('tip')
    return version

setup(
    license = 'Apache License 2.0 (http://www.apache.org/licenses/LICENSE-2.0)',
    long_description = long_description,
    use_vcs_version = {'version_handler': version_handler},
    zip_safe = True,

    author = 'HGVS Contributors',
    author_email = 'reecehart+hgvs@gmail.com',
    description = """HGVS Parser, Formatter, Mapper, Validator""",
    name = "hgvs",
    package_data = {'hgvs': ['_data/*']},
    packages = find_packages(),
    url = 'https://bitbucket.org/hgvs/hgvs',

    classifiers = [
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python",
        "Topic :: Database :: Front-Ends",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
    ],

    keywords = [
        'bioinformatics',
        'computational biology',
        'genome variants',
        'genome variation',
        'genomic variants',
        'genomic variation',
        'genomics',
        'hgvs',
    ],

    install_requires = [
        'psycopg2',
        'biopython',
        'bioutils',
        'multifastadb',
        'parsley',
        'recordtype',
    ],

    setup_requires = [
        'hgtools',
        'nose',
        # 'nose-timer', causes errors when installing via pip; cause not investigated
        'sphinx',
        'sphinxcontrib-fulltoc',
    ],

    tests_require = [
        'coverage',
        'unicodecsv',
    ],
)

## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/hgvs/hgvs)
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
