#!/bin/bash

set -e
set -u
set -o pipefail

  ## Upstream/downstream sequences of k bp length are locally aligned to locate TSDs
  ## by default, 6 bp are surveyed and up to 1 mismatch is allowed

  ## load lib
  source ${bin_path}/pfunlib.sh

  ## get sequence for TSD analysis
  ## pad_1, was deprecated as a user defined argument
  pad_1=0
  awk -v pad_1=${pad_1} -v pad_2=${pad_2} 'BEGIN{OFS="\t"}$2 > 0{print $1, $2-pad_2, $2+pad_1, $4"#"$1"_"$2"_"$3"_5p", "NA", $6}' LTR_blast_single_hit > tmp_01
  awk -v pad_1=${pad_1} -v pad_2=${pad_2} 'BEGIN{OFS="\t"}$2 > 0{print $1, $3-pad_1, $3+pad_2, $4"#"$1"_"$2"_"$3"_3p", "NA", $6}' LTR_blast_single_hit >> tmp_01

  grep -f hom_search tmp_01 \
    | sed 's/minus/\-/g' \
    | sed 's/plus/+/g' \
    | sort -k4,4 > LTR_blast_single_hit

  bedtools_getfasta sample.fasta-hm LTR_blast_single_hit LTR_blast_single_hit.fa

  ## n of LTRs / batch
  ## which is set by -r
  fasta_to_one LTR_blast_single_hit.fa > LTR_blast_single_hit.fa_fasta_to_one
  n_lines=$(wc -l < LTR_blast_single_hit.fa_fasta_to_one \
              | awk -v batch_size=${batch_size} '{if(int($0/batch_size) < $0/batch_size){print int($0/batch_size)+1}else{print int($0/batch_size)}}' \
              | awk '{if($0 % 2 != 0){print $0+1}else{print}}')
  split -l ${n_lines} LTR_blast_single_hit.fa_fasta_to_one pw_alg_batch_

  while read -r pw_alg_batch_i_; do
    tr ' ' '\n' < $pw_alg_batch_i_ > tmp_$pw_alg_batch_i_
    mv tmp_$pw_alg_batch_i_ $pw_alg_batch_i_
  done < <(ls pw_alg_batch_*)

  find $PWD -type f -name 'pw_alg_batch_*' | xargs -P ${batch_size} -n 1 ${bin_path}/mining_03_aux_01.sh

  cat tmp.*/ltr_tsd > ltr_tsd

  test_rm_ith "pw_alg_batch_*"
  test_rm_ith "tmp.*/*"
  rmdir tmp.*/

  grep -c "Smith-Waterman local alignment of sequences" water.log | awk '{print $0/2 " hits aligned"}'
  (grep -v -c "Smith-Waterman local alignment of sequences" water.log || true) | awk '$0 > 0{print "WARNING: water err!"}'


