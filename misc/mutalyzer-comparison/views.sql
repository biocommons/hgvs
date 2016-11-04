-- filter out bogus records
-- \? filters variants that have question marks in them
-- [ACGT][ACGT]> filters multinucleotide substitutions (which don't parse)
-- \d> filters records that don't have a reference nucleotide (insertion?)

create or replace view reece.snp_hgvs_filtered as
select * from reece.snp_hgvs
where hgvs_name !~ '\?' and hgvs_name !~ '[ACGT][ACGT]>' and hgvs_name !~ '\d>'
;

create or replace view reece.acmgmr_dbsnp_v as
select T.hgnc,SH.*
from reece.gene_set GS
join uta_20140210.transcript T on GS.hgnc=T.hgnc
join reece.snp_hgvs_filtered SH on T.ac=split_part(SH.hgvs_name,':',1)
where set_id=1 and T.origin_id=3 ;

create or replace view reece.acmgmr_dbsnp_tests_v as
select C.hgnc,C.snp_id,G.hgvs_name as g_hgvs,C.hgvs_name as c_hgvs
from reece.acmgmr_dbsnp_v C
join reece.snp_hgvs_filtered G on C.snp_id=G.snp_id and G.hgvs_name ~ '^NC_0000' and split_part(G.hgvs_name,':',1) in (
     'NC_000001.10','NC_000002.11','NC_000003.11','NC_000004.11','NC_000005.9',
     'NC_000006.11','NC_000007.13','NC_000008.10','NC_000009.11','NC_000010.10',
     'NC_000011.9','NC_000012.11','NC_000013.10','NC_000014.8','NC_000015.9',
     'NC_000016.9','NC_000017.10','NC_000018.9','NC_000019.9','NC_000020.10',
     'NC_000021.8','NC_000022.10','NC_000023.10','NC_000024.9')
;

create materialized view reece.acmgmr_dbsnp_tests_mv as select * from reece.acmgmr_dbsnp_tests_v;
