.SUFFIXES:
.PHONY: FORCE
.DELETE_ON_ERROR:

SHELL:=/bin/bash -e -o pipefail


default:
	@echo "Ain't no stinkin' $@ target" 1>&2; exit 1


############################################################################
## fetch "bermuda" data for transcript alignment issues
bermuda.tsv:
	psql -P footer=off -A -F$$'\t' -d uta \
	-c "select hgnc,tx_ac,s_status,b_status,sb_se_i_eq from bermuda.bermuda_data_mv where alt_ac~'^NC_0000'" \
	>$@.tmp
	mv $@.tmp $@


############################################################################
## build ACMG must report test set
SNP_HGVS.sql.gz:
	# download $@ from
	# https://bioq.saclab.net/query/downloads.php?db=bioq_dbsnp_human_137
	# and put it in $@
	@exit 1
# load snp_hgvs table
snp_hgvs-load.log: snp_hgvs.sql SNP_HGVS.sql.gz
	(cat $<; gzip -cdq $(word 2,$^) | bin/snp_hgvs-to-sql | bin/chunk) \
	| time psql -d uta 2>&1 \
	| tee $@.tmp
	mv $@.tmp $@
views.log: views.sql
	psql -d uta -f $< >$@.tmp
	mv $@.tmp $@
acmgmr.tsv: views.log
	psql -d uta -P footer=off -A -F$$'\t' -c 'select * from reece.acmgmr_dbsnp_tests_mv order by 1,3' >$@.tmp
	mv $@.tmp $@
acmgmr-g: acmgmr.tsv.gz
	gzip -cdq <$< | tail -n+2 | cut -f3 | sort -u >$@
acmgmr-c: acmgmr.tsv.gz
	gzip -cdq <$< | tail -n+2 <$< | cut -f4 | sort -u >$@


############################################################################
## run and compare g-to-c 
%-hgvs.tsv: %
	bin/hgvs-g-to-c -f $< >$@.tmp 2>$@.log
	mv $@.tmp $@

%-mzr.tsv: %
	# manually upload $< to
	# https://mutalyzer.nl/batchPositionConverter and save result
	# as $@
	@exit 1

%.cmp: %-hgvs.tsv %-mzr.tsv bermuda.tsv.gz
	./bin/hgvs-mutalyzer-compare -vv -eid -H $(word 1,$^) -M $(word 2,$^) -B $(word 3,$^) >$@.tmp 2>$@.log
	mv $@.tmp $@


############################################################################
## Utilities

%.gz:: %
	gzip -cq <$< >$@.tmp
	mv $@.tmp $@


############################################################################
## cleanup

.PHONY: clean cleaner cleanest
clean:
cleaner: clean
	rm *.mm *.mm.ac
cleanest: cleaner
	rm *.out
	rm minors.ac NLxdi.ac


############################################################################
## SCRAPS

.PHONY: garcia gcp
garcia gcp: %: %.out %.mm %.mm.ac


NLxdi.ac:
	psql -At -duta -c "select tx_ac from bermuda.bermuda_data_mv where sb_se_i_eq and s_status='NLxdi';" | sort -u >$@.tmp
	mv $@.tmp $@

minors.ac:
	psql -At -duta -c "select tx_ac from bermuda.bermuda_data_mv where sb_se_i_eq and s_minor;" | sort -u >$@.tmp
	mv $@.tmp $@

gcp.in:
	perl -lne 'next if m/^#/; print $$& while m/N[CM]\S+/g'  ../../tests/data/gcp/*.tsv | sort -u >$@.tmp
	mv $@.tmp $@

gcp-%.in: gcp.in
	perl -ne 'print if rand() < 0.$*' <$< >$@.tmp
	mv $@.tmp $@

garcia.in:
	perl -lne 'next if m/^#/; print $$& while m/N[CM]\S+/g'  ../../tests/data/garcia.tsv | sort -u >$@.tmp
	mv $@.tmp $@



%.out: %.in
	time ./hgvs-mutalyzer-compare -v -f $< 2>&1 | tee $@.tmp
	mv $@.tmp $@

%-NLxdi.out: %.in NLxdi.ac
	time ./hgvs-mutalyzer-compare -v -f $< -a $(word 2,$^) 2>&1 | tee $@.tmp
	mv $@.tmp $@

%-minors.out: %.in minors.ac
	time ./hgvs-mutalyzer-compare -v -f $< -a $(word 2,$^) 2>&1 | tee $@.tmp
	mv $@.tmp $@

%.mm: %.out
	egrep -B1 'mismatches' <$< >$@.tmp
	mv $@.tmp $@

%.mm.ac: %.mm
	perl -lne 'if (m/mismatches/) { print $$& while m/NM_[^:]+/g }' <$< | sort -u >$@.tmp
	mv $@.tmp $@
	perl -n0e 'print(join(",",map {"\x27$$_\x27"} split) )' <$@; echo


