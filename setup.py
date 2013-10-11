from setuptools import setup

with open('doc/description.rst') as f:
    long_description = f.read()


setup(
    name = "hgvs",
    use_hg_version = True,

    packages = ['hgvs'],
    package_data={'hgvs': ['grammar.txt']},

    url = 'https://bitbucket.org/invitae/hgvs',
    description = """HGVS Parser""",
    long_description = long_description,

    author = 'InVitae Keyboard Monkeys',
    author_email = 'reece+hgvs@invitae.com',
    license = 'MIT',

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
