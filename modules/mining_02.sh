#!/bin/bash

set -e
set -u
set -o pipefail

  ## To identify truncated TEs, both flanking sequences of each putative soloLTRs are compared to both internal domain termini of all intact elements provided
  ## The length of the region surveyed is 200 bp by default

  ## load lib
  source ${bin_path}/pfunlib.sh

  ## extract flanking sequences of $qlenght bp
  awk -v qlength=${qlength} 'BEGIN{OFS="\t"}{if($0 ~ "plus"){print $1, $2-qlength, $2, $4"#"$1"_"$2"_"$3"_5p", "NA", "+"}else{print $1, $2-qlength, $2, $4"#"$1"_"$2"_"$3"_5p", "NA", "-"}}' LTR_blast_single_hit \
    | cat - <(awk -v qlength=${qlength} 'BEGIN{OFS="\t"}{if($0 ~ "plus"){print $1, $3, $3+qlength, $4"#"$1"_"$2"_"$3"_3p", "NA", "+"}else{print$1, $3, $3+qlength, $4"#"$1"_"$2"_"$3"_3p", "NA", "-"}}' LTR_blast_single_hit) \
    | awk 'BEGIN{OFS="\t"}{if($2 > 0){print $0}}' > LTR_blast_single_hit_int_dom.bed

  bedtools_getfasta sample.fasta-hm LTR_blast_single_hit_int_dom.bed LTR_blast_single_hit_int_dom.fa


  ## BLAST search
  test_rm_ith hom_search

  ## if TEsorter_cls == "T"; then homology search is carried out by lineages.
  awk -F/ '$0 ~ ">"{print $2}' sample.intact.fa--split | sort | uniq | while read -r cls; do

    grep -P -A1 "/${cls}/" sample.intact.fa--split | grep -v '\-\-' > sample.intact.fa--split_${cls/\//_}
    makeblastdb -in sample.intact.fa--split_${cls/\//_} -dbtype nucl >> blast.log

    blastn -task 'blastn' -query <(grep -A1 '_5p' LTR_blast_single_hit_int_dom.fa | grep -P -A1 "/${cls}" | grep -v '\-\-') -db sample.intact.fa--split_${cls/\//_} -outfmt "6 qseqid qlen sseqid sstart send length pident qstart qend bitscore evalue" > blast_hit_5.fa-blast &
    pid_1=$!

    blastn -task 'blastn' -query <(grep -A1 '_3p' LTR_blast_single_hit_int_dom.fa | grep -P -A1 "/${cls}" | grep -v '\-\-') -db sample.intact.fa--split_${cls/\//_} -outfmt "6 qseqid qlen sseqid sstart send length pident qstart qend bitscore evalue" > blast_hit_3.fa-blast &
    pid_2=$!

    wait $pid_1
    wait $pid_2

    overlap="0.80"
    identity="80"
    awk -v identity="${identity}" -v overlap="${overlap}" 'BEGIN{OFS="\t"}{if(sqrt(($9 - $8)^2) >  $2*overlap && $7 > identity){print $0, "T"}}' <(cat blast_hit_5.fa-blast blast_hit_3.fa-blast) \
      | sort -g -k2,7r \
      | sed 's/_.p//g' \
      | awk '!seen[$1]++' >> hom_search

    test_rm_ith "sample.intact.fa--split_${cls/\//_}"
    test_rm_ith "blast_hit_*.fa-blast"

  done

  ## keep PASS list
  awk '$0 ~ ">" && $0 ~ "_5p"{print substr($0,2,length($0)-4)}' LTR_blast_single_hit_int_dom.fa \
    | grep -v -f <(awk '{print $1}' hom_search) - > tmp_01
  mv tmp_01 hom_search

  ## number of hits kept after hom assessment
  wc -l < hom_search | tally_c F "homology" - >> hit_count.txt


  ## add data into LTR_BLAST_overlap.txt file
  awk '{print $1" PASS"}' hom_search \
    | sort -k 1,1 \
    | join -1 4 -2 1 -a 1 -a 2 <(awk 'BEGIN{OFS="\t"}{$4=$4"#"$1"_"$2-1"_"$3; print}' LTR_BLAST_overlap.txt | sort -k4,4) - \
    | awk 'BEGIN{OFS="\t"}{if(NF > 8){print $0}else{print $0" FAIL"}}' | tr ' ' '\t' > tmp_01
  mv tmp_01 LTR_BLAST_overlap.txt
