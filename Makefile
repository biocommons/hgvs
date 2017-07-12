# Makefile for Python project

.DELETE_ON_ERROR:
.PHONY: FORCE
.PRECIOUS :
.SUFFIXES:

SHELL:=/bin/bash -o pipefail
SELF:=$(firstword $(MAKEFILE_LIST))

PKG=hgvs
PKGD=$(subst .,/,${PKG})
TEST_DIRS:=tests
DOC_TESTS:=doc hgvs ./README.rst


############################################################################
#= BASIC USAGE
default: help

#=> help -- display this help message
help:
	@sbin/makefile-extract-documentation "${SELF}"


############################################################################
#= SETUP, INSTALLATION, PACKAGING

#=> venv: make a Python 3 virtual environment
.PHONY: venv
venv:
	#pyvenv venv
	virtualenv venv
	source venv/bin/activate; \
	pip install --upgrade pip setuptools

#=> setup: setup/upgrade packages *in current environment*
.PHONY: setup
setup: etc/develop.reqs etc/install.reqs
	if [ -s $(word 1,$^) ]; then pip install --upgrade -r $(word 1,$^); fi
	if [ -s $(word 2,$^) ]; then pip install --upgrade -r $(word 2,$^); fi

#=> devready: create venv, install prerequisites, install pkg in develop mode
.PHONY: devready
devready:
	make venv && source venv/bin/activate && make setup develop
	@echo '#############################################################################'
	@echo '###  Do not forget to `source venv/bin/activate` to use this environment  ###'
	@echo '#############################################################################'


#=> develop: install package in develop mode
#=> install: install package
#=> bdist bdist_egg bdist_wheel build sdist: distribution options
.PHONY: bdist bdist_egg bdist_wheel build build_sphinx sdist install develop
develop:
	pip install -e .
bdist bdist_egg bdist_wheel build sdist install: %:
	python setup.py $@




############################################################################
#= TESTING
# see test configuration in setup.cfg


#=> test: execute tests
.PHONY: test
test:
	python setup.py pytest --addopts="--cov=hgvs -m 'not extra' ${TEST_DIRS}"

#=> test-docs: execute tests
.PHONY: test-docs
test-docs:
	python setup.py pytest --addopts="--cov=hgvs -m 'not extra' ${DOC_TESTS}"


#=> test-* -- run tests with specified tag
test-%:
	python setup.py pytest --addopts="--cov=hgvs -m ${*} ${TEST_DIRS}"

#=> tox -- run all tox tests
tox:
	tox

tox-%:
	tox -- -m $* $(TEST_DIRS)


############################################################################
#= UTILITY TARGETS

# N.B. Although code is stored in github, I use hg and hg-git on the command line
#=> reformat: reformat code with yapf and commit
.PHONY: reformat
reformat:
	@if hg sum | grep -qL '^commit:.*modified'; then echo "Repository not clean" 1>&2; exit 1; fi
	@if hg sum | grep -qL ' applied'; then echo "Repository has applied patches" 1>&2; exit 1; fi
	yapf -i -r "${PKGD}" tests
	hg commit -m "reformatted with yapf"

#=> docs -- make sphinx docs
.PHONY: docs
docs: develop
	# RTD makes json. Build here to ensure that it works.
	make -C doc html json

############################################################################
#= CLEANUP

#=> clean: remove temporary and backup files
.PHONY: clean
clean:
	find . \( -name \*~ -o -name \*.bak \) -print0 | xargs -0r rm
	-make -C doc $@
	-make -C examples $@

#=> cleaner: remove files and directories that are easily rebuilt
.PHONY: cleaner
cleaner: clean
	rm -fr .cache *.egg-info build dist doc/_build htmlcov
	find . \( -name \*.pyc -o -name \*.orig -o -name \*.rej \) -print0 | xargs -0r rm
	find . -name __pycache__ -print0 | xargs -0r rm -fr
	/bin/rm -f examples/.ipynb_checkpoints
	/bin/rm -f hgvs.{dot,svg,png,sfood}
	-make -C doc $@
	-make -C examples $@

#=> cleanest: remove files and directories that require more time/network fetches to rebuild
.PHONY: cleanest
cleanest: cleaner
	rm -fr .eggs .tox venv
	-make -C doc $@
	-make -C examples $@

#=> pristine distclean: above, and delete anything unknown to mercurial
pristine distclean: cleanest
	if [ -d .hg ]; then hg st -inu0 | xargs -0r /bin/rm -fv; fi


############################################################################
#= HGVS specific

changelog:
	make -C doc/changelog 0.4.rst

#=> changelog
doc/source/changelog.rst: CHANGELOG
	./sbin/clog-txt-to-rst <$< >$@

#=> code-check -- check code with flake8
code-check:
	flake8 hgvs test --output-file=$@.txt

#=> docs-aux -- make generated docs for sphinx
docs-aux:
	make -C misc/railroad doc-install
	make -C examples doc-install

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


## <LICENSE>
## Copyright 2016 Source Code Committers
## 
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
## 
##     http://www.apache.org/licenses/LICENSE-2.0
## 
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
