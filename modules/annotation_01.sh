#!/bin/bash

set -e
set -u
set -o pipefail

  ## all soloLTR descriptors are gathered and final gff3 file is generated

  ## load lib
  source ${bin_path}/pfunlib.sh

  ## keep hits whose TSD <= n mismatch
  tsd_count ltr_tsd ${mismatch} > tmp_01
  mv tmp_01 ltr_tsd

  cut -f2 ltr_tsd \
    | awk -F '[#]' 'BEGIN{OFS="\t"}{split($3,arr,"_"); print arr[1]"_"arr[2], arr[3], arr[4], $3"#"$2}' \
    | sort -k 1,1 -k 2,2n > ltr_tsd.bed


  ## get soloLTR seq
  bedtools_getfasta sample.fasta ltr_tsd.bed sample.intact.fa--soloLTR

  ## n of LTRs / batch
  ## which is set by -r
  n_lines=$(wc -l < ltr_tsd.bed \
              | awk -v batch_size=${batch_size} '{if(int($0/batch_size) < $0/batch_size){print int($0/batch_size)+1}else{print int($0/batch_size)}}' \
              | awk '{if($0 % 2 != 0){print $0+1}else{print}}')
  split -l ${n_lines} ltr_tsd.bed anno_batch_

  find $PWD -type f -name 'anno_batch_*' | xargs -P ${batch_size} -n 1 ${bin_path}/annotation_01_aux_01.sh
  cat tmp.*/soloLTR.gff3 > soloLTR.gff3

  hits_=$(wc -l soloLTR.gff3 | awk '{print $1}')
  sort -k1,1 -k4,4n soloLTR.gff3 \
    | awk -v hits_=$hits_ 'BEGIN{OFS="\t"}{if(length(NR) < length(hits_)){leading_ = length(hits_)-length(NR); id_ = sprintf("%0*d", leading_, 0)NR; sub("Name=soloLTR_.*;Classification","Name=soloLTR_"id_";Classification",$NF); print}else{sub("Name=soloLTR_.*;Classification","Name=soloLTR_"NR";Classification",$NF); print}}' > tmp000_
  mv tmp000_ soloLTR.gff3

  test_rm_ith "anno_batch_*"
  test_rm_ith "tmp.*/*"
  rmdir tmp.*/

  ## extract soloLTR sequences with proper header
  awk 'BEGIN{OFS="\t"}{split($9,arr,"=|;"); print $1, $4-1, $5, arr[2], $6, $7}' soloLTR.gff3 > soloLTR.bed
  bedtools_getfasta sample.fasta soloLTR.bed soloLTR.fa
