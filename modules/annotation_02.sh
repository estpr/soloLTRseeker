#!/bin/bash

set -e
set -u
set -o pipefail


  ## final output files are sorted and working directory cleaned

  ## load lib
  source ${bin_path}/pfunlib.sh

  ## add TSD data into LTR_BLAST_overlap.txt file
  grep -f <(awk '$0 !~ /^##/{print $1"_"$4-1"_"$5}' soloLTR.gff3) ltr_tsd \
    | awk 'BEGIN{OFS="\t"}{print $2"\tPASS"}' \
    | sort -k1,1 \
    | join -1 1 -2 1 -a 1 -a 2 <(sort -k1,1 LTR_BLAST_overlap.txt | tr ' ' '\t') - \
    | awk '{if(NF > 9){print}else{print $0" FAIL"}}' | tr ' ' '\t' > tmp_01
  mv tmp_01 LTR_BLAST_overlap.txt

  ## add TSD seq
  grep -f <(awk '$0 !~ /^##/{print $1"_"$4-1"_"$5}' soloLTR.gff3) ltr_tsd \
    | awk 'BEGIN{OFS="\t"}{print $2, $3"_"$4}' \
    | sort -k1,1 \
    | join -1 1 -2 1 -a 1 -a 2 <(sort -k1,1 LTR_BLAST_overlap.txt | tr ' ' '\t') - \
    | awk '{if(NF > 10){print $1, $2, $3, $4, $5, $6, $7, $NF, $8, $9, $10}else{print $1, $2, $3, $4, $5, $6, $7, "NA", $8, $9, $10}}' | tr ' ' '\t' > tmp_01
  mv tmp_01 LTR_BLAST_overlap.txt

  ## ann locus type
  ## complete type refers to those loci whose all hits are covered by the longest one
  ## conversely, partial tag describes those loci for which both ends bear overhanging nucleotides in reference to the whole overlapping region
  sort -k 7,7 -k 3,3n -k 4,4nr LTR_BLAST_overlap.txt \
    | awk 'BEGIN{prev_locus = $7; prev_end = $4}{if(prev_locus == $7 && prev_end < $4){print $0,"partial";prev_end = $4}else if(prev_locus == $7 && prev_end >= $4){print $0, "complete";prev_end = prev_end}else{print $0, "complete"; prev_locus = $7; prev_end = $4}}' \
    | tr -s ' ' \
    | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $NF, $8, $9, $10, $11}' | tr ' ' '\t' > tmp_01
  mv tmp_01 LTR_BLAST_overlap.txt

  cut -f7,8 LTR_BLAST_overlap.txt | sort -rd -k 2,2 | awk '!a[$1]++' > locus.list
  map_val locus.list LTR_BLAST_overlap.txt 7 8 > tmp_01
  mv tmp_01 LTR_BLAST_overlap.txt

  ## number of hits +TSD
  wc -l soloLTR.gff3 | tally_c T "TSD" - >> hit_count.txt

  ## number of hits +TSD_no_mismatch
  grep -c -P "_0$" soloLTR.gff3 | tally_c F "TSD_no_mismatch" - >> hit_count.txt

  ## number of hits +TSD_in_partial_locus
  (awk '$NF == "PASS" && $NF == $(NF-1) && $NF == $(NF-2) && $(NF-3) == "tpartial"{print $7}' LTR_BLAST_overlap.txt || true) \
    | sort \
    | uniq \
    | wc -l | tally_c F "TSD_in_partial_locus" - >> hit_count.txt


