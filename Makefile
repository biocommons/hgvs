# Makefile for Python project

.DELETE_ON_ERROR:
.PHONY: FORCE
.PRECIOUS:
.SUFFIXES:

ifeq ("","$(shell command -v zsh)")
$(error "zsh not found; you must install zsh first")
endif

SHELL:=zsh -eu -o pipefail -o null_glob
SELF:=$(firstword $(MAKEFILE_LIST))

VE_DIR=venv

SRC_DIRS:=src
TEST_DIRS:=tests
DOC_TESTS:=docs hgvs ./README.rst

# TESTING sources
# UTA_DB_URL must be accessible at test time (for now)
export HGVS_CACHE_MODE=
_UTAPW=anonymous
export UTA_DB_URL=postgresql://anonymous:${_UTAPW}@uta.biocommons.org:5432/uta/uta_20241220
export HGVS_SEQREPO_URL=http://localhost:5000/seqrepo


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
${VE_DIR}:
	python3 --version
	python3 -mvenv $@; \
	source $@/bin/activate; \
	python3 -m ensurepip --upgrade; \
	pip install --upgrade pip setuptools uv wheel

#=> develop: install package in develop mode
.PHONY: develop
develop:
	uv pip install -e ".[dev]"
	pre-commit install

#=> install: install package
.PHONY: install
install:
	uv pip install "."

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
.PHONY: test test-code test-docs
test:
	HGVS_CACHE_MODE=run pytest --cov ${SRC_DIRS}
test-code:
	pytest ${TEST_DIRS}
test-docs:
	pytest ${DOC_TESTS}
stest:
	pytest -vvv -s -k ${t}
test-%:
	pytest -m '$*' ${TEST_DIRS}

#=> test-learn: add new data to test data cache
test-learn:
	HGVS_CACHE_MODE=learn pytest -s
#=> test-relearn: destroy and rebuild test data cache
test-relearn:
	rm -fr tests/data/cache-py3.hdp tests/cassettes
	HGVS_CACHE_MODE=learn pytest -s
#=> test-relearn-iteratively: destroy and rebuild test data cache (biocommons/hgvs#760)
# temporary workaround for https://github.com/biocommons/hgvs/issues/760
test-relearn-iteratively:
	rm -fr tests/data/cache-py3.hdp tests/cassettes
	find tests/ -name 'test*.py' | HGVS_CACHE_MODE=learn xargs -tn1 -- pytest --no-cov -x -s

#=> tox -- run all tox tests
tox:
	tox

#=> cqa: execute code quality tests
cqa:
	ruff format --check
	ruff check

############################################################################
#= UTILITY TARGETS

#=> reformat: reformat code
.PHONY: reformat
reformat:
	ruff check --fix
	ruff format

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
	-make -C docs $@
	-make -C examples $@

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
