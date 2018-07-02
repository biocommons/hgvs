This directory contains examples of challenges when projecting
variants in the vicinity of alignment discrepancies.

discrepencies := substitution, insertion, and deletion alignment differences

projected variant locations may be within, exactly equal to, partially
overlap, or subsume the discrepancy

See `./run-tests` for defined tests. An excerpt:

example:
```
(3.6) snafu$ ./run-tests
(3.6) snafu$ tail +2 output/1.1.3.dev12+nc0fa49f31de4.d20180617 | column -t -s$'\t'
disc_type  loc_type  var_n                             expected                                     got (if failed)
sub        within    NM_183425.2:n.41G>C               NC_000020.11:g.57391438=
sub        within    NM_183425.2:n.41G>T               NC_000020.11:g.57391438C>T
gdel       within    NM_007121.5:n.797A>T              NC_000019.10:g.50378563_50378564insATC       NC_000019.10:g.50378564_50378563delinsT
...
```
