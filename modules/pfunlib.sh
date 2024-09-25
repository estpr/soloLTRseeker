#!/bin/bash

set -e
set -u
set -o pipefail

  ## extract RC sequences for downstream analysis
  bedtools_getfasta () {
    if [ -f ${1}.fai ];then
      rm -I ${1}.fai
    fi
    printf  "${1}.fai has been removed\n \n" >> test_rm_ith.log 
    bedtools getfasta -name -s -fi ${1} -bed ${2} -fo ${3}
    sed -i 's/:.*//g' ${3}
    sed -i 's/(.)//g' ${3}
  }


  ## find_locus function is intended to group BLAST hits located in the *same* locus under the same identifier
  ## locus is defined here as a feature whose length is extended as long as the next adjacent hit overlaps the previous one, forming an uninterrupted unit
  find_locus () {
    awk -v overlap=${overlap} 'BEGIN{OFS="\t"}{if(sqrt(($9 - $8)^2) >  $2*overlap){print $0, "T"}}' ${1} \
      | awk 'BEGIN {OFS="\t"}{if($4 > $5){print $1, $2, $3, $5, $4, $6, $7, $8, $9, $10, $11, $12, $13, $14}else{print $0}}' \
      | sort -k 3,3 -k 4,4n -k 5,5nr \
      | awk 'BEGIN {prev_end = $5; chr = $3; id = "locus_"NR}{if(prev_end >= $4 && chr == $3 && prev_end >= $5){print $0"\tNA"; prev_end = prev_end; chr = $3; id = "locus_"NR}else if(prev_end >= $4 && chr == $3 && prev_end < $5){print $0"\tNA"; prev_end = $5; chr = $3; id = "locus_"NR}else{print $0"\t"id; prev_end = $5; chr = $3; id = "locus_"NR}}' \
      | awk 'BEGIN {locus = "locus_0"} $13 ~ "NA" {$13 = locus} {locus = $13} 1' \
      | awk 'BEGIN {OFS="\t"}{print $3, $4, $5, $1, $10, $7, $11, $5-$4, $13}'
  }


  ## parse_alg extracts relevant information from pw-alignments saved in markx2 and markx3 format (${seq})
  parse_alg () {
    seq=${2}
    grep -v '#' ${1} \
      | grep -v '^$\|^ \|>' \
      | tr '\n' '\t' \
      | awk -F '[\t]' -v seq=${seq} '{if(seq==1){split($1,a," ");print a[2]}else{split($2,b," ");print b[length(b)]}}'
  }


  ## replace adapted chr names by true ones
  remap_chr () {
    chr=${1}
    key=${2}
    sed -i "s/${key}\([^0-9]\)/${chr}\1/g" ${3}
  }


  ## change tmp file names
  rename_output () {
    mv ${1} ${sample}_${1}
  }


  ## reverse complement dna seq
  reverse_complement () {
    grep -o '.' ${1} \
      | awk 'NR == FNR{a[$1]=$2;next} {$1=a[$1]}1' <(printf "A\tT\nT\tA\nC\tG\nG\tC\n") - \
      | tac - \
      | tr -d '\n'
  }


  ## annotate # of hits at every step of the analysis
  tally_c () {
    awk -v sample="${sample}" -v run="${run}" -v module="${module}" -v header="${1}" -v step="${2}" 'BEGIN{OFS="\t"}{if(header == "T"){print "## "module": \n"sample, run, step, $1}else{print sample, run, step, $1}}' ${3}
  }


  ## remove files in a controlled way
  test_rm_ith () {
    counter=0
    for file in ${1};do
      if [ -f ${file} ];then
        counter=$[$counter +1]
      fi
    done

    ## test how many files will be treated
    if [[ $counter -lt 20 ]];then
      printf "$counter files will be rm:\n" >> test_rm_ith.log
    else
      printf "Too many args. Aborted!" >> test_rm_ith.log
      exit 1
      ## return 1 might be used too as long as you capture it;i.e., (test_rm_ith) || exit $?
    fi

    for file in ${1};do
      if [ -f ${file} ];then
        printf "file $file rm\n" >> test_rm_ith.log
        rm -I ${file}
      fi
    done
    printf " \n" >> test_rm_ith.log
  }


  ## to keep only TSD that match the threshods 
  tsd_count () {
    awk -v mismatch=${2} '{if($5 >= 4 && $5 <= 6 && $5 == $6 && $8 <= mismatch){print}}' ${1} \
      | awk '{if($5 == 4 && $8 != 0){ next }{print}}' \
      | awk -v mismatch=${mismatch} '{if($5 == 4 && mismatch == 0){ next }{print}}'
  }


  ## count unique entries in a tab
  var_count () {
    sort ${1} | uniq -c | awk '{$1=$1};1' | sort -k1,1nr | tr ' ' '\t'
  }

