.SUFFIXES :
.PRECIOUS :
.PHONY : FORCE
.DELETE_ON_ERROR:

SHELL:=/bin/bash -o pipefail
SELF:=$(firstword $(MAKEFILE_LIST))

ifdef LOCAL_UTA
export UTA_DB_URL=postgresql://localhost/uta
endif

############################################################################
#= BASIC USAGE
default: help

#=> help -- display this help message
help:
	@sbin/extract-makefile-documentation "${SELF}"

############################################################################
#= INSTALLATION/SETUP

#=> ve -- create a virtualenv
VE_DIR:=ve
VE_MAJOR:=1
VE_MINOR:=10
VE_PY_DIR:=virtualenv-${VE_MAJOR}.${VE_MINOR}
VE_PY:=${VE_PY_DIR}/virtualenv.py
${VE_PY}:
	curl -sO  https://pypi.python.org/packages/source/v/virtualenv/virtualenv-${VE_MAJOR}.${VE_MINOR}.tar.gz
	tar -xvzf virtualenv-${VE_MAJOR}.${VE_MINOR}.tar.gz
	rm -f virtualenv-${VE_MAJOR}.${VE_MINOR}.tar.gz
${VE_DIR}: ${VE_PY} 
	${SYSTEM_PYTHON} $< ${VE_DIR} 2>&1 | tee "$@.err"
	/bin/mv "$@.err" "$@"


############################################################################
#= UTILITY TARGETS

#=> develop, build_sphinx, sdist, upload_sphinx
bdist bdist_egg build build_sphinx develop install sdist upload_docs: %:
	python setup.py $*

# sphinx docs needs to be able to import packages
build_sphinx: develop

#=> docs -- make sphinx docs
docs: # build_sphinx

test-with-coverage: test-setup
	python setup.py nosetests --with-xunit --with-coverage --cover-erase --cover-package=hgvs --cover-html --exclude test_nightly

test-all-with-coverage: test-setup
	python setup.py nosetests --with-xunit --with-coverage --cover-erase --cover-package=hgvs --cover-html 


#=> upload
upload:
	python setup.py bdist_egg sdist upload


#=> lint -- run lint, flake, etc
# TBD


#=> test-setup -- prepare to run tests
test-setup:

#=> test -- run tests
test: install
	python setup.py nosetests --with-xunit

#=> continuous integration tests -- target for jenkins (and now travis, drone, or other providers)
ci-test jenkins:
	make ve \
	&& source ve/bin/activate \
	&& make install \
	&& make test-with-coverage \
	&& make docs

ci-test-nightly jenkins-nightly:
	make cleanest \
	&& make ve \
	&& source ve/bin/activate \
	&& make install \
	&& make test-all-with-coverage \
	&& make docs


############################################################################
hgvs/data/seguid-acs.json.gz:
	gzip -cdq human.protein.faa.gz* | ./sbin/fasta-seguid | gzip -cq >$@


############################################################################
#= CLEANUP
.PHONY: clean cleaner cleanest pristine
#=> clean: clean up editor backups, etc.
clean:
	find . -name \*~ -print0 | xargs -0r /bin/rm
#=> cleaner: above, and remove generated files
cleaner: clean
	find . -name \*.pyc -print0 | xargs -0r /bin/rm -f
	/bin/rm -fr distribute-* *.egg *.egg-info *.tar.gz nosetests.xml
	-make -C doc clean
#=> cleanest: above, and remove the virtualenv, .orig, and .bak files
cleanest: cleaner
	find . \( -name \*.orig -o -name \*.bak \) -print0 | xargs -0r /bin/rm -v
	/bin/rm -fr build bdist dist sdist ve virtualenv*
#=> pristine: above, and delete anything unknown to mercurial
pristine: cleanest
	hg st -un0 | xargs -0r echo /bin/rm -fv

## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/invitae/hgvs)
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
