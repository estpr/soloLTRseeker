#!/bin/bash

set -e
set -u
set -o pipefail


  ## directory is cleaned and final output files renamed

  ## load lib
  source ${bin_path}/pfunlib.sh


  ## map chromosome names
  cut -f2 chr_map.txt | sort | uniq | while read -r map; do
    chr=$(awk -v map=${map} '$2 == map{print $1}' chr_map.txt || true)
    key=$(awk -v map=${map} '$2 == map{print $2}' chr_map.txt || true)

    remap_chr ${chr} ${key} LTRRT_lengths.txt || true
    remap_chr ${chr} ${key} LTR_BLAST_overlap.txt || true
    remap_chr ${chr} ${key} soloLTR.gff3 || true
    remap_chr ${chr} ${key} soloLTR.fa || true
  done

  ## append headers
  awk 'BEGIN{OFS="\t"}{split($1,arr,"#"); print $2, $3, $4, arr[1]"#"arr[2], $5, $6, $7, $8, $9, $10, $11, $12}' LTR_BLAST_overlap.txt \
    | cat <(echo "## chr start end queryID strand piden locus locus_type TSD_alg BLAST_single_hit homology TSD") - > tmp_01
  mv tmp_01 LTR_BLAST_overlap.txt

  echo "## queryID FL_length lLTR_length rLTR_length intdom_length motif length_diff tg_ca cd_hit" \
    | cat -  <(cat LTRRT_lengths.txt) > tmp_01
  mv tmp_01 LTRRT_lengths.txt


  ## rename outputfiles
  if [ -z ${target_chr+x} ]; then
    export target_chr="wga"
  fi

  rename_output hit_count.txt "${target_chr}"
  rename_output LTRRT_lengths.txt "${target_chr}"
  rename_output LTR_BLAST_overlap.txt "${target_chr}"
  rename_output soloLTR.gff3 "${target_chr}"
  rename_output soloLTR.fa "${target_chr}"

  ## cleaning directory

  test_rm_ith "sample.fasta*"
  test_rm_ith "LTR*"
  test_rm_ith "ltr_*"
  test_rm_ith "tmp_01*"
  test_rm_ith "hom_search"
  test_rm_ith "soloLTR.bed"
  test_rm_ith "locus.list"

  cut -f1 chr_map.txt | while read -r chr; do
    test_rm_ith "genome_dir/${chr}.fasta"
  done
  rmdir genome_dir

  test_rm_ith "chr_map.txt"

  if [[ "${TEsorter_cls}" == "T" ]]; then
    awk 'NR>1{print $2}' sample.cls | awk '!a[$0]++' | while read -r cls; do
      test_rm_ith "sample.intact.fa--split_${cls}*"
    done
  fi

  test_rm_ith "sample.cls"
  test_rm_ith "sample.intact*"
