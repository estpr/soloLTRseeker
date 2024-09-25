#!/bin/bash

set -e
set -u
set -o pipefail

  ## if -t argument is switched on, otherwise this step is skipped

  ## load lib
  source ${bin_path}/pfunlib.sh

  ## generate fasta file for LTRRTs
  awk 'BEGIN{OFS="\t"}$11 == "fl_LTRRT"{$3=$10;print $1, $4-1, $5, $3, $6, $7}' sample.intact.gff3 \
    | grep -f <(grep -P "\tPASS" LTRRT_lengths.txt | cut -f1) - > sample.intact.bed

  bedtools_getfasta sample.fasta sample.intact.bed sample.intact.fa

  TEsorter sample.intact.fa -db rexdb-plant -p 10 -nolib > sample.intact.fa--TEsorter.log 2>&1

  ## keep hits with consistent domain annotation
  grep -v 'none$' sample.intact.fa.rexdb-plant.cls.tsv \
    | paste -d '\t' - <(grep -v 'none$' sample.intact.fa.rexdb-plant.cls.tsv \
                         | cut -f7 \
                         | awk -F '[ |]' '{for(i=2;i<=NF;i+=2){printf $i" "}; printf "\n" }' \
                         | awk '{ for (i=2; i<=NF; ++i) if ($i != $1) { print "F"; next } print "T" }') \
    | awk -F '[\t]' 'BEGIN{OFS="\t"}{if($8 == "T"){print $1,$3"_"$4}}' > sample.cls

  ## add cls data into sample.intact.gff3 file
  ## unclassified LTRRTs will be assigned a generic NA
  sort -k 1,1 sample.cls \
    | join -1 1 -2 10 -a 1 -a 2 - <(sort -k 10,10 sample.intact.gff3) \
    | awk 'BEGIN{OFS="\t"}{if(NF == 12){print $3, $4, $5, $6, $7, $8, $9, $10, $11, $1"/"$2, $12}else{print $2, $3, $4, $5, $6, $7, $8, $9, $10, $1"/NA", $11}}' \
    | tr ' ' '\t' > tmp_01
  mv tmp_01 sample.intact.gff3

  ## add cls data into LTRRT_lengths.txt file
  sort -k 1,1 sample.cls \
    | join -1 1 -2 1 -a 1 -a 2 <(sort -k 1,1 LTRRT_lengths.txt) - \
    | awk 'BEGIN{OFS="\t"}{if(NF > 6){$1=$1"/"$NF;print}else{$1=$1"/NA";print}}' \
    | cut -f1-6 \
    | grep -v '^#' > tmp_01
  mv tmp_01 LTRRT_lengths.txt

  ## remove TEsorter output
  test_rm_ith "sample.intact.fa.rexdb-plant*"

  printf "\n\nTEsorter classification done.\n\n"
