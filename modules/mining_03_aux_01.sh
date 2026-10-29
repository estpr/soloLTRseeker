#!/bin/bash

set -e
set -u
set -o pipefail

  ## Upstream/downstream sequences of k bp length are locally aligned to locate TSDs
  ## by default, 6 bp are surveyed and up to 1 mismatch is allowed
  ## parrallel run

  ## load lib
  source ${bin_path}/pfunlib.sh

  pw_alg_batch_=${1}

  tmpdir_=$(mktemp -d tmp.XXXXXXXXXX)
  outdir_=$PWD

  ## align seq
  test_rm_ith "ltr_tsd"

  unset ltr_tsd_arr
  grep '>' $pw_alg_batch_ | sed 's/>\|_.p//g' | sort | uniq | while read -r ltr; do

    grep -A1 "${ltr}" LTR_blast_single_hit.fa | grep -v '\-\-' > $tmpdir_/ltr.fa

    water <(head -n2 $tmpdir_/ltr.fa) <(tail -n2 $tmpdir_/ltr.fa) -aformat markx2 -outfile $tmpdir_/ltr.fa.alg.i -gapopen 16 -gapextend 4 2>> water.log
    water <(head -n2 $tmpdir_/ltr.fa) <(tail -n2 $tmpdir_/ltr.fa) -aformat markx3 -outfile $tmpdir_/ltr.fa.alg.ii -gapopen 16 -gapextend 4 2>> water.log

    sample="MusMusc_GCA_056825265.1.fna.F2B"
    ltr_tsd_arr[0]=$(echo ${sample})
    ltr_tsd_arr[1]=$(echo ${ltr})
    ltr_tsd_arr[2]=$(parse_alg $tmpdir_/ltr.fa.alg.i 1)
    ltr_tsd_arr[3]=$(parse_alg $tmpdir_/ltr.fa.alg.ii 2)
    ltr_tsd_arr[4]=$(echo ${#ltr_tsd_arr[2]})
    ltr_tsd_arr[5]=$(echo ${#ltr_tsd_arr[3]})
    ltr_tsd_arr[6]=$(parse_alg $tmpdir_/ltr.fa.alg.i 2 | awk 'BEGIN{OFS="\t"}{print gsub(/\./,"",$1)}')
    ltr_tsd_arr[7]=$(( ${ltr_tsd_arr[4]} - ${ltr_tsd_arr[6]} ))
    ltr_tsd_arr[8]=$(parse_alg $tmpdir_/ltr.fa.alg.i 2 | awk -F '[^.]' -vORS=';' '{i=1; while(i<=NF) {print length($i); i++} } END {printf "\n"}')

    echo ${ltr_tsd_arr[@]} | tr ' ' '\t' >> $tmpdir_/ltr_tsd

    unset ltr_tsd_arr
    test_rm_ith $tmpdir_/ltr.fa.alg.i
    test_rm_ith $tmpdir_/ltr.fa.alg.ii
    test_rm_ith $tmpdir_/ltr.fa
  done
