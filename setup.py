import os
import sys

from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

# full path appears to be required for old (0.6.x?) versions of setuptools
root_dir = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(root_dir, 'doc/source/description.rst')) as f:
    long_description = f.read()

setup(
    author = 'InVitae Keyboard Monkeys',
    license = 'Apache License 2.0 (http://www.apache.org/licenses/LICENSE-2.0)',
    long_description = long_description,
    use_hg_version = True,
    zip_safe = False,

    author_email = 'reece+hgvs@invitae.com',
    description = """HGVS Parser""",
    name = "hgvs",
    package_data = {'hgvs': ['data/*']},
    packages = find_packages(),
    url = 'https://bitbucket.org/invitae/hgvs',

    install_requires = [
        'bdi',
        'biopython',
        'parsley',
        'recordtype',
        'uta',
        ],

    setup_requires = [
        'hgtools',
        'nose',
        'sphinx',
        'sphinx_rtd_theme',
        'sphinxcontrib-httpdomain',
        ],

    tests_require = [
        'coverage',
        'nose',
        ],
    )

## <LICENSE>
## Copyright 2014 UTA Contributors (https://bitbucket.org/invitae/uta)
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
