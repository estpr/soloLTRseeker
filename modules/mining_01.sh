#!/bin/bash

set -e
set -u
set -o pipefail

  ## blastn run is used to identify candidate loci across the target genome
  ## A minimum of 0.99 query coverage is defined by default.

  ## load lib
  source ${bin_path}/pfunlib.sh

  ## mask reference genome with gff3 input file
  bedtools maskfasta -fi sample.fasta -bed <(grep 'fl_LTRRT' sample.intact.gff3) -fo sample.fasta-hm

  ## BLAST
  makeblastdb -in sample.fasta-hm -dbtype nucl >> blast.log

  ## n of LTRs / batch
  n_lines=$(wc -l < sample.intact.fa--LTR \
              | awk -v batch_size=${batch_size} '{if(int($0/batch_size) < $0/batch_size){print int($0/batch_size)+1}else{print int($0/batch_size)}}' \
              | awk '{if($0 % 2 != 0){print $0+1}else{print}}')
  split -l ${n_lines} sample.intact.fa--LTR batch_

  if [[ ${parallel} == "T" ]]; then

    echo "blastn"
    pids=()
    while read -r file_batch; do
      blastn -task 'blastn' -query batch_${file_batch} -db sample.fasta-hm -evalue 1e-6 -qcov_hsp_perc ${overlap##*.} -perc_identity 80 -outfmt "6 qseqid qlen sseqid sstart send length pident qstart qend sstrand evalue" > blast_batch_${file_batch} &
      pids+=($!)
    done < <(ls -htl batch_* | awk -F '[_]' '{print $NF}')

    status=()
    for pid in "${pids[@]}"; do
      wait "${pid}"
      status+=($?)
    done

    echo " "
    for i in "${!status[@]}"; do
      echo "job $i exited with ${status[$i]}"
    done

    cat blast_batch_* > LTR_blast

  else

    blastn -task 'blastn' -query sample.intact.fa--LTR -db sample.fasta-hm -evalue 1e-6 -qcov_hsp_perc ${overlap##*.} -perc_identity 80 -num_threads ${batch_size} -outfmt "6 qseqid qlen sseqid sstart send length pident qstart qend sstrand evalue" > LTR_blast

  fi


  ## remove batch files
  test_rm_ith "blast_batch_*"
  test_rm_ith "batch_*"

  cp LTR_blast check_blast
  ## create LTR_BLAST_overlap.txt file
  find_locus LTR_blast \
    | sort -k 1,1 -k 2,2n -k 3,3n \
    | tr ' ' '\t' \
    | cut -f1-6,9 > LTR_BLAST_overlap.txt


  ## number of hits overlap > ${overlap}
  wc -l LTR_BLAST_overlap.txt | tally_c T "BLAST_overlap" - >> hit_count.txt


  ## keep single hits with greater p.iden
  find_locus LTR_blast \
    | sort -k 8,8nr -k 6,6nr \
    | awk '!a[$9]++' \
    | sort -k 1,1 -k 2,2n -k 3,3n \
    | awk 'BEGIN{OFS="\t"}{print $1, $2-1, $3, $4, "NA", $5}' > LTR_blast_single_hit


  ## add data into LTR_BLAST_overlap.txt file
  awk 'BEGIN{OFS="\t"}{print $1"&"$2+1"&"$3"&"$4, $5, $6, $7, $9}' LTR_blast_single_hit \
   | sort -k1,1 \
   | join -1 1 -2 1 -a 1 -a 2 <(awk 'BEGIN{OFS="\t"}{print $1"&"$2"&"$3"&"$4, $5, $6, $7, $9}' LTR_BLAST_overlap.txt | sort -k1,1) - \
   | awk '{if(NF > 5){print $1, $2, $3, $4, "PASS"}else{print $1, $2, $3, $4, "FAIL"}}' \
   | awk '!a[$1]++' \
   | sed 's/ \|&/\t/g' > tmp_01
  mv tmp_01 LTR_BLAST_overlap.txt


  ## number of hits after collapsing each locus hits to single representatives
  wc -l LTR_blast_single_hit | tally_c F "BLAST_single_hit" - >> hit_count.txt

