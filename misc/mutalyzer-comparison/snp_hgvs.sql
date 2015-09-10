drop table reece.snp_hgvs;

create table reece.snp_hgvs (
  snp_id integer,
  hgvs_name text not null,
  source text default null,
  upd_time timestamp default null
);

create index snp_hgvs_snp_id on reece.snp_hgvs(snp_id);

