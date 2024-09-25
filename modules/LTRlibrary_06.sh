#!/bin/bash

set -e
set -u
set -o pipefail

  ## Both lLTR sides are quereyed for TG..CA pattern conservation
  ## By default, 1 mismatch is allowed

  ## load lib
  source ${bin_path}/pfunlib.sh

  ## check TG..CA pattern
  awk 'BEGIN{OFS = "\t"}{if($0 ~ ">"){ORS = "\t"; print substr($0,2,length($0))}else{ORS = "\n"; up_dnt = substr($0,1,2); down_dnt = substr($0,length($0)-1,2); print toupper(up_dnt), toupper(down_dnt)}}' sample.intact.fa--LTR > ltr_motif

  grep -A1 -f <(grep $'.G\tCA\|T.\tCA\|TG\t.A\|TG\tC.' ltr_motif | cut -f1) sample.intact.fa--LTR | grep -v '\-\-' > tmp_01
  mv tmp_01 sample.intact.fa--LTR


  ## number of lLTR kept after TG..CA
  grep -c '>' sample.intact.fa--LTR | tally_c F "TG..CA" - >> hit_count.txt


  ## add cls data into LTRRT_lengths.txt file
  awk -F/ '$0 ~ ">"{print substr($0,2,length($0))" PASS"}' sample.intact.fa--LTR \
    | sort -k1,1 \
    | join -1 1 -2 1 -a 1 -a 2 <(sort -k1,1 LTRRT_lengths.txt) - \
    | awk 'BEGIN{OFS="\t"}{if(NF==6){print $0" FAIL"}else{print}}' > tmp_01
  mv tmp_01 LTRRT_lengths.txt

  awk 'BEGIN{OFS="\t"}{print $1, $2"_"$3}' ltr_motif \
    | sort -k1,1 \
    | join -1 1 -2 1 -a 1 -a 2 <(sort -k1,1 LTRRT_lengths.txt) - \
    | awk 'BEGIN{OFS="\t"}{if(NF==7){print $1, $2, $3, $4, $5, "NA", $6, $7}else{print $1, $2, $3, $4, $5, $8, $6, $7}}' > tmp_01
  mv tmp_01 LTRRT_lengths.txt
