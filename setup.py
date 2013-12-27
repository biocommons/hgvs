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
        'nose',
        'sphinx',
        'sphinx-pypi-upload',
        'sphinx_rtd_theme',
        'sphinxcontrib-httpdomain',

        'biopython',
        'parsley',
        'recordtype',

        # Non-PyPI dependencies -- requires dependency_links below
        # N.B. These do not work via pip install, but do work with python setup.py install
        'bdi',
        'uta',
        ],

    dependency_links = [
        # for non-PyPI dependencies
        'hg+ssh://hg@bitbucket.org/locusdevelopment/bdi#egg=bdi',
        'hg+ssh://hg@bitbucket.org/locusdevelopment/uta#egg=uta',
        ],

    setup_requires = [
        'coverage',
        'hgtools',
        ],

    )
