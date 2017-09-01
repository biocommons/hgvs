from setuptools import setup, find_packages
from pip.req import parse_requirements

with open('doc/description.txt') as f:
    long_description = f.read()

install_reqs = parse_requirements('requirements.txt')
# reqs is a list of requirement
# e.g. ['django==1.5.1', 'mezzanine==1.4.6']
reqs = [str(ir.req) for ir in install_reqs]

setup(
    license = 'Apache License 2.0 (http://www.apache.org/licenses/LICENSE-2.0)',
    long_description = long_description,
    use_scm_version = True,
    zip_safe = True,

    author = 'HGVS Contributors',
    author_email = 'hgvs-discuss@gmail.com',
    description = """HGVS Parser, Formatter, Mapper, Validator""",
    name = "hgvs",
    package_data = {'hgvs': ['_data/*']},
    packages = find_packages(),
    url = 'https://bitbucket.org/biocommons/hgvs',

    entry_points={
        'console_scripts': [
            'hgvs-shell = hgvs.shell:shell'
            ],
        },

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


    setup_requires = [
        'setuptools_scm',
        'nose',
        'sphinx==1.1.3',
        'sphinx_rtd_theme',
        'sphinxcontrib-fulltoc>=1.1',
        'wheel',
    ],

    install_requires=reqs,

    tests_require = [
        'coverage',
        'unicodecsv',
    ],

)

## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
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
