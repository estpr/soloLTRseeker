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
    | sed 's/plus/+/g' > LTR_blast_single_hit

  bedtools_getfasta sample.fasta-hm LTR_blast_single_hit LTR_blast_single_hit.fa

  ## align seq
  test_rm_ith "ltr_tsd"

  unset ltr_tsd_arr
  grep '>' LTR_blast_single_hit.fa | sed 's/>\|_.p//g' | sort | uniq | while read -r ltr; do

    grep -A1 "${ltr}" LTR_blast_single_hit.fa | grep -v '\-\-' > ltr.fa

    water <(head -n2 ltr.fa) <(tail -n2 ltr.fa) -aformat markx2 -outfile ltr.fa.alg.i -gapopen 16 -gapextend 4 2>> water.log
    water <(head -n2 ltr.fa) <(tail -n2 ltr.fa) -aformat markx3 -outfile ltr.fa.alg.ii -gapopen 16 -gapextend 4 2>> water.log

    ltr_tsd_arr[0]=$(echo ${sample})
    ltr_tsd_arr[1]=$(echo ${ltr})
    ltr_tsd_arr[2]=$(parse_alg ltr.fa.alg.i 1)
    ltr_tsd_arr[3]=$(parse_alg ltr.fa.alg.ii 2)
    ltr_tsd_arr[4]=$(echo ${#ltr_tsd_arr[2]})
    ltr_tsd_arr[5]=$(echo ${#ltr_tsd_arr[3]})
    ltr_tsd_arr[6]=$(parse_alg ltr.fa.alg.i 2 | awk 'BEGIN{OFS="\t"}{print gsub(/\./,"",$1)}')
    ltr_tsd_arr[7]=$(( ${ltr_tsd_arr[4]} - ${ltr_tsd_arr[6]} ))
    ltr_tsd_arr[8]=$(parse_alg ltr.fa.alg.i 2 | awk -F '[^.]' -vORS=';' '{i=1; while(i<=NF) {print length($i); i++} } END {printf "\n"}')

    echo ${ltr_tsd_arr[@]} | tr ' ' '\t' >> ltr_tsd

    unset ltr_tsd_arr
    test_rm_ith ltr.fa.alg.i
    test_rm_ith ltr.fa.alg.ii
    test_rm_ith ltr.fa

  done

  grep -c "Smith-Waterman local alignment of sequences" water.log | awk '{print $0/2 " hits aligned"}'
  (grep -v -c "Smith-Waterman local alignment of sequences" water.log || true) | awk '$0 > 0{print "WARNING: water err!"}'


