#import os, sys
#root_dir = os.path.dirname(__file__)
#sys.path[0:0] = [os.path.join(root_dir ,'lib','python')]

from distribute_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

# fetch __version__
with open('hgvs/version.py') as f:
    exec(f.read())

with open('doc/description.rst') as f:
    long_description = f.read()

#pkg_dir = 'lib/python'
setup(
    author = 'InVitae Keyboard Monkeys',
    author_email='reece.hart@invitae.com',
    description = """HGVS Parser and Formatter""",
    license = 'MIT',
    long_description = long_description,
    name = "hgvs",
    #package_dir = {'': pkg_dir},
    #packages = find_packages(pkg_dir),
    url = 'https://bitbucket.org/invitae/hgvs',
    version = __version__,
    zip_safe = True,
    install_requires = [
        'nose',
        'sphinx',
        'sphinx-pypi-upload',
        ],    
)
