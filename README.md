# *hgvs* - manipulate biological sequence variants according to Human Genome Variation Society recommendations

The *hgvs* package provides a Python library to parse, format, validate,
normalize, and map sequence variants according to [Variation
Nomenclature](http://varnomen.hgvs.org/) (aka Human Genome Variation
Society) recommendations.

Specifically, the hgvs package focuses on the subset of the HGVS
recommendations that precisely describe sequence-level variation
relevant to the application of high-throughput sequencing to clinical
diagnostics. The package does not attempt to cover the full scope of
HGVS recommendations. Please refer to
[issues](https://github.com/biocommons/hgvs/issues) for limitations.

## Information

[![rtd](https://img.shields.io/badge/docs-readthedocs-green.svg)](http://hgvs.readthedocs.io/) [![changelog](https://img.shields.io/badge/docs-changelog-green.svg)](https://hgvs.readthedocs.io/en/stable/changelog/)  [![getting_help](https://img.shields.io/badge/!-help%20me-red.svg)](https://hgvs.readthedocs.io/en/stable/getting_help.html)  [![GitHub license](https://img.shields.io/github/license/biocommons/hgvs.svg)](https://github.com/biocommons/hgvs/blob/main/LICENSE)  [![binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/biocommons/hgvs/main?filepath=examples)

## Latest Release

[![GitHub tag](https://img.shields.io/github/tag/biocommons/hgvs.svg)](https://github.com/biocommons/hgvs) [![pypi_rel](https://img.shields.io/pypi/v/hgvs.svg)](https://pypi.org/project/hgvs/)

## Development

[![coveralls](https://img.shields.io/coveralls/github/biocommons/hgvs.svg)](https://coveralls.io/github/biocommons/hgvs) [![issues](https://img.shields.io/github/issues-raw/biocommons/hgvs.svg)](https://github.com/biocommons/hgvs/issues)
[![GitHub Open Pull Requests](https://img.shields.io/github/issues-pr/biocommons/hgvs.svg)](https://github.com/biocommons/hgvs/pull/) [![GitHub license](https://img.shields.io/github/contributors/biocommons/hgvs.svg)](https://github.com/biocommons/hgvs/graphs/contributors/) [![GitHub stars](https://img.shields.io/github/stars/biocommons/hgvs.svg?style=social&label=Stars)](https://github.com/biocommons/hgvs/stargazers) [![GitHub forks](https://img.shields.io/github/forks/biocommons/hgvs.svg?style=social&label=Forks)](https://github.com/biocommons/hgvs/network)

## Features

- Parsing is based on formal grammar.
- An easy-to-use object model that represents most variant types
  (SNVs, indels, dups, inversions, etc) and concepts (intronic
  offsets, uncertain positions, intervals)
- A variant normalizer that rewrites variants in canonical forms and
  substitutes reference sequences (if reference and transcript
  sequences differ)
- Formatters that generate HGVS strings from internal representations
- Tools to map variants between genome, transcript, and protein
  sequences
- Reliable handling of regions genome-transcript discrepancies
- Pluggable data providers support alternative sources of transcript
  mapping data
- Extensive automated tests, including those for all variant types and
  \"problematic\" transcripts
- Easily installed using remote data sources. Installation with local
  data sources is straightforward and completely obviates network
  access

## Citation

Wang M, Callenberg KM, Dalgleish R, Fedtsov A, Fox N, Freeman PJ, et al.
<br/>**hgvs: A Python package for manipulating sequence variants using HGVS nomenclature: 2018 Update.**
<br/>Hum Mutat. 2018. [doi:10.1002/humu.23615](https://doi.org/10.1002/humu.23615)

## Important Notes

- **You are encouraged to** [browse
  issues](https://github.com/biocommons/hgvs/issues). All known issues
  are listed there. Please report any issues you find.
- **Use a pip package specification to stay within minor releases.**
  For example, `hgvs>=1.5,<1.6`. hgvs uses [Semantic
  Versioning](http://semver.org/).

## Installing HGVS Locally

**Important:** For more detailed installation and configuration
instructions, see the [HGVS readthedocs](https://hgvs.readthedocs.io/)

### Prerequisites

    libpq
    python3
    postgresql

Examples for installation:

MacOS :

    brew install libpq
    brew install python3
    brew install postgresql@14

Ubuntu :

    sudo apt install gcc libpq-dev python3-dev

### Installation Steps

By default, hgvs uses remote data sources, which makes
installation easy. If you would like to use local instances of the data sources, see the [readthedocs](https://hgvs.readthedocs.io/).

1. Create a virtual environment using your preferred method.

    **Example:**

        python3 -m venv venv

2. Run the following commands in your virtual environment:

        source venv/bin/activate
        pip install --upgrade setuptools
        pip install hgvs

See [Installation
instructions](http://hgvs.readthedocs.org/en/stable/installation.html)
for details, including instructions for installing [Universal Transcript
Archive (UTA)](https://github.com/biocommons/uta/) and
[SeqRepo](https://github.com/biocommons/biocommons.seqrepo/) locally.

## Examples and Usage

See [examples](https://github.com/biocommons/hgvs/tree/main/examples) and [readthedocs](https://hgvs.readthedocs.io/) for usage.

## Contributing

The hgvs package is a community effort. Please see
[Contributing](http://hgvs.readthedocs.org/en/stable/contributing.html) to get
started in submitting source code, tests, or documentation. Thanks for getting
involved!

### Testing

Existing tests use a cache that is committed with the repo to ensure that tests do not require external networking.  To develop new tests, which requires loading the cache, you should install UTA and Seqrepo (and the rest service) locally.

    docker compose --project-name biocommons -f $PWD/misc/docker-compose.yml up

IMPORTANT: Loading the test caches is currently hampered b
[#551](https://github.com/biocommons/hgvs/issues/551),
[#760](https://github.com/biocommons/hgvs/issues/760), and
[#761](https://github.com/biocommons/hgvs/issues/761). To load reliably, use
`make test-relearn-iteratively` for now.

## See Also

Other packages that manipulate HGVS variants:

- [pyhgvs](https://github.com/counsyl/hgvs)
- [Mutalyzer](https://mutalyzer.nl/)
