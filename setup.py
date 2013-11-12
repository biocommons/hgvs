import os
import sys

from setuptools import setup

# full path appears to be required for old (0.6.x?) versions of setuptools
root_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(root_dir, 'doc/description.rst')) as f:
    long_description = f.read()

setup(
    author = 'InVitae Keyboard Monkeys',
    author_email = 'reece+hgvs@invitae.com',
    description = """HGVS Parser""",
    license = 'MIT',
    long_description = long_description,
    name = "hgvs",
    package_data={'hgvs': ['grammar.txt']},
    packages = ['hgvs'],
    url = 'https://bitbucket.org/invitae/hgvs',
    use_hg_version = True,
    zip_safe = True,

    install_requires = [
        'nose',
        'parsley',
        'recordtype',
        'sphinx',
        'sphinx-pypi-upload',
        ],

    setup_requires = [
        'hgtools',
        ]
)
