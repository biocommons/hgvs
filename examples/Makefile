.PHONY: FORCE
.SUFFIXES:
.DELETE_ON_ERROR:

nbs=$(subst .ipynb,,$(wildcard *.ipynb))

default all: $(addsuffix .rst,${nbs})

%.rst: %.ipynb
	ipython nbconvert --to rst $<
