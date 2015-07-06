#!/bin/sh

# Demonstration of a bug with samtools that results in incorrect sequences
# being returned when the source fasta has been bgzip'd.

# WARNING: This script Creates and removes files in current directory
# invoke in /tmp or somewhere safe

# For samtools 0.1.19-96b5f2294a on x86_64 GNU/Linux (Ubuntu 13.10), running
# script returns:
#### snafu$ ~/projects/hgvs/misc/samtools-bug.sh 
#### + samtools
#### + head -3
#### 
#### Program: samtools (Tools for alignments in the SAM format)
#### Version: 0.1.19-96b5f2294a
#### + hexdump -C lf.fa
#### 00000000  3e 73 31 0a 41 43 47 54  0a 3e 73 32 0a 54 47 43  |>s1.ACGT.>s2.TGC|
#### 00000010  41 0a                                             |A.|
#### 00000012
#### + hexdump -C crlf.fa
#### 00000000  3e 73 31 0d 0a 41 43 47  54 0d 0a 3e 73 32 0d 0a  |>s1..ACGT..>s2..|
#### 00000010  54 47 43 41 0d 0a                                 |TGCA..|
#### 00000016
#### + samtools faidx lf.fa s2
#### >s2
#### TGCA
#### + samtools faidx lf.fa.gz s2
#### >s2
#### >s1A
#### + samtools faidx crlf.fa s2
#### >s2
#### TGCA
#### + samtools faidx crlf.fa.gz s2
#### >s2
#### >s1A

## The lf, crlf tests demonstrate that the bug is attributable to
## compression and not linefeed flavor.


############################################################################

# create the base test file
cat <<EOF >lf.fa
>s1
ACGT
>s2
TGCA
EOF

# unix to dos newline conversion
perl -p0e 's/\n/\r\n/g' <lf.fa >crlf.fa

bgzip <lf.fa >lf.fa.gz
bgzip <crlf.fa >crlf.fa.gz

# sanity check: compress-decompress recovers original
bgzip -d <lf.fa.gz | cmp lf.fa -
bgzip -d <crlf.fa.gz | cmp crlf.fa -

samtools faidx lf.fa
samtools faidx lf.fa.gz
samtools faidx crlf.fa
samtools faidx crlf.fa.gz


set -x
samtools 2>&1 | head -3

hexdump -C lf.fa
hexdump -C crlf.fa

samtools faidx lf.fa s2
samtools faidx lf.fa.gz s2
samtools faidx crlf.fa s2
samtools faidx crlf.fa.gz s2

