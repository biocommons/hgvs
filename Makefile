# Makefile for Python project

.DELETE_ON_ERROR:
.PHONY: FORCE
.PRECIOUS:
.SUFFIXES:

SHELL:=/bin/bash -e -o pipefail
SELF:=$(firstword $(MAKEFILE_LIST))

PY_VERSION:=$(shell python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
VE_DIR=venv/${PY_VERSION}

TEST_DIRS:=tests
DOC_TESTS:=docs hgvs ./README.rst


############################################################################
#= BASIC USAGE
default: help

#=> help -- display this help message
help:
	@sbin/makefile-extract-documentation "${SELF}"


############################################################################
#= SETUP, INSTALLATION, PACKAGING

#=> devready: create venv, install prerequisites, install pkg in develop mode
.PHONY: devready
devready:
	make ${VE_DIR} && source ${VE_DIR}/bin/activate && make develop
	@echo '#################################################################################'
	@echo '###  Do not forget to `source ${VE_DIR}/bin/activate` to use this environment  ###'
	@echo '#################################################################################'

#=> venv: make a Python 3 virtual environment
venv: ${VEDIR}
${VE_DIR}: venv/%:
	python$* -mvenv $@; \
	source $@/bin/activate; \
	python -m ensurepip --upgrade; \
	pip install --upgrade pip setuptools wheel

#=> develop: install package in develop mode
.PHONY: develop install
develop:
	@if [ -z "$${VIRTUAL_ENV}" ]; then echo "Not in a virtual environment; see README.md" 1>&2; exit 1; fi
	pip install -e .[dev]

#=> install: install package
install:
	pip install .

#=> build: make sdist and wheel
.PHONY: build
build: %:
	python -m build


############################################################################
#= TESTING
# see test configuration in setup.cfg

#=> test: execute tests
#=> test-code: test code (including embedded doctests)
#=> test-docs: test example code in docs
#=> stest: test a specific test or set of tests, also useful using ipdb
#=> test-/tag/ -- run tests marked with /tag/
# TODO: rationalize tags
# find tests -name \*.py | xargs perl -ln0e 'while (m/@pytest.mark.(\w+)/g) {print $1 if not $seen{$1}++}'  | sort
# => extra fx issues mapping models normalization parametrize pnd quick regression validation
.PHONY: test test-code test-docs stest
test:
	pytest
test-code:
	pytest ${TEST_DIRS}
test-docs:
	pytest ${DOC_TESTS}
stest:
	pytest -vvv -s -k ${t}
test-%:
	pytest -m '$*' ${TEST_DIRS}

#=> test-learn: build test data cache
#=> test-relearn: destroy and rebuild test data cache
test-learn:
	# The appropriate environment should be set up using misc/docker-compose.yml
	UTA_DB_URL=postgresql://anonymous@localhost:5432/uta/uta_20180821 \
	HGVS_SEQREPO_URL=http://localhost:5000/seqrepo \
	HGVS_CACHE_MODE=learn \
	pytest -s
test-relearn:
	rm -f tests/data/cache-py3.hdp
	make test-learn

#=> tox -- run all tox tests
tox:
	tox


############################################################################
#= UTILITY TARGETS

#=> reformat: reformat code with yapf and commit
.PHONY: reformat
reformat:
	@if ! git diff --cached --exit-code >/dev/null; then echo "Repository not clean" 1>&2; exit 1; fi
	black src tests
	isort src tests
	git commit -a -m "reformatted with black and isort"

#=> docs -- make sphinx docs
.PHONY: docs
docs: develop
	# RTD makes json. Build here to ensure that it works.
	make -C docs html json

############################################################################
#= CLEANUP

#=> clean: remove temporary and backup files
.PHONY: clean
clean:
	rm -frv **/*~ **/*.bak
	-make -C docs $@
	-make -C examples $@

#=> cleaner: remove files and directories that are easily rebuilt
.PHONY: cleaner
cleaner: clean
	rm -frv .cache build dist docs/_build
	rm -frv **/__pycache__
	rm -frv **/*.egg-info
	rm -frv **/*.pyc
	rm -frv **/*.orig
	rm -frv **/*.rej

#=> cleanest: remove files and directories that are more expensive to rebuild
.PHONY: cleanest
cleanest: cleaner
	rm -frv .eggs .tox venv
	-make -C docs $@
	-make -C examples $@

#=> distclean: remove untracked files and other detritus
.PHONY: distclean
distclean: cleanest
	git clean -df


############################################################################
#= HGVS specific

#=> code-check -- check code with flake8
code-check:
	flake8 hgvs test --output-file=$@.txt

#=> docs-aux -- make generated docs for sphinx
docs-aux:
	make -C misc/railroad docs-install
	make -C examples docs-install

#=> hgvs.svg -- import graph; requires snakefood
hgvs.sfood:
	sfood hgvs >"$@.tmp"
	/bin/mv "$@.tmp" "$@"
hgvs.dot: hgvs.sfood
	sfood-graph -p <$< >"$@.tmp"
	/bin/mv "$@.tmp" "$@"
hgvs.svg: hgvs.dot
	dot -Tsvg <$< >"$@.tmp"
	/bin/mv "$@.tmp" "$@"


# <LICENSE>
# Copyright 2018 HGVS Contributors (https://github.com/biocommons/hgvs)
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
