This directory contains examples of challenges when projecting
variants in or around regions of alignment discrepancies.  Some are
bugs in hgvs and others represent intrinsic biological difficulties
with some alignments.

There are two basic alignment discrepancies we are concerned with:
conservative sequence changes and indels. 


                      variant type
                   sub           ins          del

how variant location relates to alignment discrepancy
within				   
exact
partial
covers

Also, we should probably consider both substitution and gap alignment discrepancies

That is, I think that there are probably 3 types ⊗ 4 location
relationships ⊗ 2 discrepancy types = 24 cases to consider.


## Gap alignment discrepancy


```
           57391435  57391442
                |      |
NC_000020.11    GCACGTCG
NM_183425.2     GCAGCTCG
                |  ^^  |
               38      45
```

```
            50378561  50378566
                |       |
NC_000019.10    GGA---AAC
NM_007121.5     GGAAACAAC
                |       |
               793     801
```
```
           149779572  149779580
                |       |
NC_000007.14    TGACAGCCC
NM_198455.2     TGA---CCC
                |       |
              1113     1118
```



## Substitution alignment discrepancy

|  type   |             variant                |  hgvs projected variant (n->g, normalize=False) |     note     |
|---------|------------------------------------|-------------------------------------------------|--------------|
| within  |  NM_183425.2:n.41G>C               |   NC_000020.11:g.57391438=                      |  reasonable  |
|         |  NM_183425.2:n.41G>T               |   NC_000020.11:g.57391438C>T                    |  reasonable  |
|         |  NM_183425.2:n.41=                 |   NC_000020.11:g.57391438C>G                    |  reasonable  |
|         |  NM_183425.2:n.41del               |   NC_000020.11:g.57391438del                    |  reasonable  |
|         |  NM_183425.2:n.41_42insT           |   NC_000020.11:g.57391438_57391439insT          |  reasonable  |
| exact   |  NM_183425.2:n.41_42delinsCG       |   NC_000020.11:g.57391438_57391439=             |  reasonable  |
|         |  NM_183425.2:n.41_42delinsTAA      |   NC_000020.11:g.57391438_57391439delinsTAA     |  reasonable  |
|         |  NM_183425.2:n.41_42del            |   NC_000020.11:g.57391438_57391439del           |  reasonable  |
| partial |  NM_183425.2:n.40_41delinsTC       |   NC_000020.11:g.57391437_57391438delinsTC      |  reasonable  |
|         |  NM_183425.2:n.40_41delinsTT       |   NC_000020.11:g.57391437_57391438delinsTT      |  reasonable  |
|         |  NM_183425.2:n.40_41del            |   NC_000020.11:g.57391437_57391438del           |  reasonable  |
|         |  NM_183425.2:n.40_41insT           |   NC_000020.11:g.57391437_57391438insT          |  reasonable  |
|         |  NM_183425.2:n.42_43insG           |   NC_000020.11:g.57391439_57391440insG          |  reasonable  |
| covers  |  NM_183425.2:n.40_43del            |   NC_000020.11:g.57391437_57391440del           |  reasonable  |
|         |  NM_183425.2:n.40_43delinsCG       |   NC_000020.11:g.57391437_57391440delinsCG      |  reasonable  |


|  type   |             variant                |  hgvs projected variant (n->g, normalize=False) |              expected variant               |     note     |
|---------|------------------------------------|-------------------------------------------------|---------------------------------------------|--------------|
| within  |  NM_007121.5:n.797A>T              |   NC_000019.10:g.50378564_50378563delinsT       |  NC_000019.10:g.50378563_50378564insATC     |  bug         |
|         |  NM_007121.5:n.797del              |   NC_000019.10:g.50378564_50378563del           |  NC_000019.10:g.50378563_50378564insAC      |  bug         |
|         |  NM_007121.5:n.796_797insT         |   NC_000019.10:g.50378564_50378563insT          |  NC_000019.10:g.50378563_50378564insATAC    |  bug         |
|         |                                    |                                                 |                                             |              |
| exact   |  NM_007121.5:n.796_798=            |   NC_000019.10:g.50378564_50378563delinsAAC     |  NC_000019.10:g.50378563_50378564insAAC     |  bug         |
|         |  NM_007121.5:n.796_798del          |   NC_000019.10:g.50378564_50378563del           |  NC_000019.10:g.50378563_50378564=          |  bug         |
|         |  NM_007121.5:n.796_798delinsTCGG   |   NC_000019.10:g.50378564_50378563delinsTCGG    |  NC_000019.10:g.50378563_50378564insTCGG    |  bug         |
|         |                                    |                                                 |                                             |              |
| partial |  NM_007121.5:n.795_796del          |   NC_000019.10:g.50378563del                    |  NC_000019.10:g.50378563delinsAC            |  bug         |
|         |  NM_007121.5:n.795_796delinsTT     |   NC_000019.10:g.50378563delinsTT               |  NC_000019.10:g.50378563delinsTTAC          |  bug         |
|         |  NM_007121.5:n.795_796insT         |   NC_000019.10:g.50378563insT                   |  NC_000019.10:g.50378563_50378564insTAAC    |  bug         |
|         |                                    |                                                 |                                             |              |
| covers  |  NM_007121.5:n.794_800del          |   NC_000019.10:g.50378562_50378565del           |                                             |  reasonable  |
|         |  NM_007121.5:n.794_800delinsTC     |   NC_000019.10:g.50378562_50378565delinsTC      |                                             |  reasonable  |


|  type   |             variant                |  hgvs projected variant (n->g, normalize=False) |              expected variant               |     note     |
|---------|------------------------------------|-------------------------------------------------|---------------------------------------------|--------------|
| within  |  NM_198455.2:n.1115_1116insT       |   NC_000007.14:g.149779574_149779578insT        |  NC_000007.14:g.149779575_149779577delinsT  |  bug         |
|         |                                    |                                                 |                                             |              |
| exact   |  NM_198455.2:n.1115_1116insCAG     |   NC_000007.14:g.149779574_149779578insCAG      |  NC_000007.14:g.149779574_149779578=        |  bug         |
|         |  NM_198455.2:n.1115_1116=          |   NC_000007.14:g.149779574_149779578delinsAC    |                                             |  reasonable  |
|         |                                    |                                                 |                                             |              |
| partial |  NM_198455.2:n.1115_1116insCA      |   NC_000007.14:g.149779574_149779578insC        |  NC_000007.14:g.149779575_149779577delinsCA |  bug         |
|         |                                    |                                                 |                                             |              |
| cover   |  NM_198455.2:n.1114_1117del        |   NC_000007.14:g.149779573_149779579del         |                                             |  reasonable  |
|         |  NM_198455.2:n.1114_1117delinsCAG  |   NC_000007.14:g.149779573_149779579delinsCAG   |                                             |  reasonable  |


