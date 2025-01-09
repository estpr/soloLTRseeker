#!/bin/bash

set -e
set -u
set -o pipefail

  ## chromosome names are adapted to avoid regex conflicts

  ## load lib
  source ${bin_path}/pfunlib.sh

  ## check viability of chromosome names
  if grep -q ">" <(awk '$0 ~ ">" && $0 ~ /[^a-zA-Z0-9>_-]/' sample.fasta) ; then
    printf "\nSpecial characters found in chromosome names. It's likely they cause problems!\n"
    printf "We'll try to adapt them removing any non-alphanumeric character...\n"

    awk '{if($0 ~ ">"){gsub(/[^a-zA-Z0-9>]/,"",$0); print}else{print}}' sample.fasta > tmp_01
    mv tmp_01 sample.fasta

    awk -F '[\t]' '{gsub(/[^a-zA-Z0-9]/,"",$1); print}' sample.intact.gff3 > tmp_01
    mv tmp_01 sample.intact.gff3
  fi

  ## generate chromosome name tab
  awk '$0 ~ ">"' sample.fasta | awk '{chr = substr($0,2,length($0)); print chr"\tchr_"NR}' > chr_map.txt


  ## split genome and rename headers
  mkdir genome_dir
  csplit -z -s -f genome_dir/chr sample.fasta "/^>/" {*}

  for file in genome_dir/chr*; do
    mv ${file} "genome_dir/$(head -1 ${file} | cut -d/ -f1 | tr -d '>').fasta"
  done

  while read -r chr; do
    key=$(grep -P "${chr}\t" chr_map.txt | cut -f2)
    awk -v key=${key} '{if(NR == 1){$0 = ">"key; print}else{print}}' genome_dir/${chr}.fasta >> temp.fasta
  done < <(cut -f1 chr_map.txt)
  mv temp.fasta sample.fasta


  ## rename chromosome in gff3 file
  map_val chr_map.txt sample.intact.gff3 1 1 \
    | awk '$3 !~ /repeat_region|target_site_duplication/' \
    | awk 'BEGIN{OFS="\t"}$0 !~ /^ /{print $0, $1"_"$4, $1"_"$5, sqrt(($5-$4)^2)}' \
    | sort -k1,1 -k4,4n -k12,12nr \
    | awk 'BEGIN{OFS="\t"}NF--' > tmp_01
  mv tmp_01 sample.intact.gff3
