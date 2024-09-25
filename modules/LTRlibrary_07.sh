#!/bin/bash

set -e
set -u
set -o pipefail

  ## LTRs are reduced to a non redundant set of representatives using cd-hit-est and a default percentage identity of 0.95
  ## Coverage cutoffs are preset to aL == 0.3 and aS == 0.9

  ## load lib
  source ${bin_path}/pfunlib.sh

  ## cluster LTR sequences
  cd-hit-est -c ${cd_hit} -aL ${aL} -aS ${aS} -g 1 -d 0 -p 1 -i sample.intact.fa--LTR -o tmp_01
  mv tmp_01 sample.intact.fa--LTR

  ## number of lLTR kept after clustering
  grep -c '>' sample.intact.fa--LTR | tally_c F "cd-hit" - >> hit_count.txt

  ## add cls data into LTRRT_lengths.txt file
  awk -F/ '$0 ~ ">"{print substr($0,2,length($0))" PASS"}' sample.intact.fa--LTR \
    | sort -k1,1 \
    | join -1 1 -2 1 -a 1 -a 2 <(sort -k1,1 LTRRT_lengths.txt) - \
    | awk 'BEGIN{OFS="\t"}{if(NF==8){print $0" FAIL"}else{print}}' | tr ' ' '\t' > tmp_01
  mv tmp_01 LTRRT_lengths.txt
