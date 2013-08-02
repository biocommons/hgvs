from distutils.core import setup

with open('hgvs/version.py') as f:
    exec(f.read())

with open('doc/description.rst') as f:
    long_description = f.read()

setup(
    name = "hgvs",
    version = __version__,
    packages = ['hgvs', 'hgvs.utils'],

    url = 'https://bitbucket.org/invitae/hgvs',
    description = """HGVS Parser and Formatter""",
    long_description = long_description,

    author = 'InVitae Keyboard Monkeys',
    author_email='reece.hart@invitae.com',
    license = 'MIT',

    zip_safe = True,

    install_requires = [
        'nose',
        'parsley',
        'recordtype',
        'sphinx',
        'sphinx-pypi-upload',
        ],    
)
