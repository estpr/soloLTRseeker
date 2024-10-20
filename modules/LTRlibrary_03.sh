#!/bin/bash

set -e
set -u
set -o pipefail


  ## min LTR length and max size divergence between counterparts are tested
  ## By default, 100 bp and 50 bp are used as thresholds, respectively.

  ## load lib
  source ${bin_path}/pfunlib.sh

  ## combine feature lengths
  awk 'BEGIN{OFS="\t"}$11=="fl_LTRRT"{print $10, $5-$4+1}' sample.intact.gff3 \
    | sort -k 1,1 \
    | paste -d '\t' - <(awk 'BEGIN{OFS="\t"}$11=="lLTR"{print $10, $5-$4+1}' sample.intact.gff3 | sort -k 1,1) <(awk 'BEGIN{OFS="\t"}$11=="rLTR"{print $10, $5-$4+1}' sample.intact.gff3 | sort -k 1,1) \
    | awk 'BEGIN{OFS="\t"}{if($1 == $3 && $1 == $5){print}else{print "unsorted file!"}}' > tmp_01

  if grep "unsorted file!" tmp_01; then
    printf "LTRRT coordinates coundn't be retrieved!"
    exit 1
  fi

  ## eval LTR_diff and min_LTR_length
  awk -v LTR_diff=${LTR_diff} -v min_LTR_length=${min_LTR_length} 'BEGIN{OFS="\t"}{if(sqrt(($4-$6)^2) > LTR_diff || $4 < min_LTR_length){print $1,$2,$4,$6,$2-($4+$6),"FAIL"}else{print $1,$2,$4,$6,$2-($4+$6),"PASS"}}' tmp_01 > LTRRT_lengths.txt

  ## number of PASS_fl_LTRRT
  grep -c 'PASS' LTRRT_lengths.txt | tally_c F "PASS_fl_LTRRT" - >> hit_count.txt

  if [[ "${TEsorter_cls}" == "F" ]]; then
    ## add NA tag into lineage cls loc
    awk 'BEGIN{OFS="\t"}{$10=$10"/NA";print}' sample.intact.gff3 > tmp_01
    mv tmp_01 sample.intact.gff3

    ## add NA tag into lineage cls loc
    awk 'BEGIN{OFS="\t"}{$1=$1"/NA";print}' LTRRT_lengths.txt > tmp_01
    mv tmp_01 LTRRT_lengths.txt
  fi
