.PHONY: FORCE
.SUFFIXES:
.DELETE_ON_ERROR:

SHELL:=/bin/bash

DOC_EXAMPLES_DIR:=../doc/source/examples

nbs=$(subst .ipynb,,$(wildcard *.ipynb))
rst_nbs:=$(addsuffix .rst,${nbs})

default all: ${rst_nbs}

%.rst: %.ipynb
	ipython nbconvert --to rst $<

doc-install: ${rst_nbs}
	rm -fr ${DOC_EXAMPLES_DIR}; mkdir ${DOC_EXAMPLES_DIR}; cp $^ ${DOC_EXAMPLES_DIR}/
