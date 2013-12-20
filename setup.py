import os
import sys

from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

# full path appears to be required for old (0.6.x?) versions of setuptools
root_dir = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(root_dir, 'doc/description.rst')) as f:
    long_description = f.read()

setup(
    author = 'InVitae Keyboard Monkeys',
    license = 'MIT',
    long_description = long_description,
    use_hg_version = True,
    zip_safe = True,

    author_email = 'reece+hgvs@invitae.com',
    description = """HGVS Parser""",
    name = "hgvs",
    package_data = {'hgvs': ['data/*']},
    packages = find_packages(),
    url = 'https://bitbucket.org/invitae/hgvs',

    install_requires = [
        'biopython',
        'nose',
        'parsley',
        'pysam',
        'recordtype',
        'sphinx',
        'sphinx-pypi-upload',
        ],

    setup_requires = [
        'hgtools',
        ],
)
