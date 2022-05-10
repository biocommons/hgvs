from setuptools import setup, find_packages
from sys import version_info

short_description = "HGVS Parser, Formatter, Mapper, Validator"
with open("docs/description.txt") as f:
    long_description = f.read()

setup(license="Apache License 2.0 (http://www.apache.org/licenses/LICENSE-2.0)",
      long_description=long_description,
      use_scm_version=True,
      zip_safe=True,
      author="HGVS Contributors",
      author_email = 'hgvs-discuss@gmail.com',
      description = short_description.replace("\n", " "),
      name="hgvs",
      package_data={"hgvs": ["_data/*"]},
      packages=find_packages(),
      url="https://github.com/biocommons/hgvs",
      entry_points={
          'console_scripts': [
              'hgvs-shell = hgvs.shell:shell'
          ],
      },
      classifiers=[
          "Development Status :: 5 - Production/Stable",
          "Intended Audience :: Developers",
          "Intended Audience :: Healthcare Industry",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: Apache Software License",
          "Operating System :: OS Independent",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python",
          "Topic :: Database :: Front-Ends",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          "Topic :: Scientific/Engineering :: Medical Science Apps.",
      ],
      keywords=[
          "bioinformatics",
          "computational biology",
          "genome variants",
          "genome variation",
          "genomic variants",
          "genomic variation",
          "genomics",
          "hgvs",
      ],
      install_requires=[
          "attrs>=17.4.0",  # https://github.com/biocommons/hgvs/issues/473
          "biocommons.seqrepo>=0.6.5",
          "bioutils>=0.4.0,<1.0",
          "configparser>=3.3.0",
          "ipython",
          "parsley",
          "psycopg2",
          "six"
      ],
      setup_requires=[
          "cython",             # required for RTD build :-)
          "pytest-runner",
          "setuptools_scm",
          "wheel",
      ],
      tests_require=[
          "pytest>=5.3",
          "pytest-cov>=2.8",
      ],
)


# <LICENSE>
# Copyright 2021 HGVS Contributors (https://github.com/biocommons/hgvs)
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
