#!/bin/bash

set -e
set -u
set -o pipefail


  ## sequences for downstream analysis are extracted with appropriate headers
  ## generate bed and fasta files

  ## load lib
  source ${bin_path}/pfunlib.sh

  ## lLTR
  awk 'BEGIN{OFS="\t"}$11=="lLTR"{print $1, $4-1, $5, $10, $6, $7}' sample.intact.gff3 \
    | grep -f <(awk '$NF == "PASS"{print $1}' LTRRT_lengths.txt) - > sample.intact.bed--LTR

  bedtools_getfasta sample.fasta sample.intact.bed--LTR sample.intact.fa--LTR


  ## internal domain ${qlength} bp edge sequences
  ## end coordinate of lLTR and start coordinate of rLTSs are used
  test_rm_ith "sample.intact.bed--split"

  awk -v qlength=${qlength} 'BEGIN{OFS="\t"}$11=="lLTR"{print $1, $5-1, $5+qlength, $10"/"$11, $6, $7}' sample.intact.gff3 \
    | grep -f <(awk '$NF == "PASS"{print $1}' LTRRT_lengths.txt) - >> sample.intact.bed--split
  awk -v qlength=${qlength} 'BEGIN{OFS="\t"}$11=="rLTR"{print $1, $4-qlength-1, $4, $10"/"$11, $6, $7}' sample.intact.gff3 \
    | grep -f <(awk '$NF == "PASS"{print $1}' LTRRT_lengths.txt) - >> sample.intact.bed--split

  bedtools_getfasta sample.fasta sample.intact.bed--split sample.intact.fa--split
