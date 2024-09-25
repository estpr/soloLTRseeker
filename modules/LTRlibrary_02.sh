#!/bin/bash

set -e
set -u
set -o pipefail


  ## gff3 input file is scanned to identify and flag fl_LTR, lLTR and rLTR

  ## load lib
  source ${bin_path}/pfunlib.sh

  ## fl_LTRRT
  awk 'NR==FNR && !/^#/{a[$10]++;b[$11]++;next}a[$10] > 1 && b[$11] > 1 && a[$10] < 3 && b[$11] < 3' sample.intact.gff3 sample.intact.gff3 \
    | awk '{print $0"\t"$1"_"$4"_"$5"#"$3"\tfl_LTRRT";}' \
    | tr ' ' '\t' \
    | cut -f1-9,12,13 > tmp_01

  ## lLTR
  awk 'NR==FNR && !/^#/{a[$10]++;b[$11]++;next}a[$10] > 1 && a[$10] < 3' sample.intact.gff3 sample.intact.gff3 \
    | awk 'NR==FNR && !/^#/{a[$10]++;b[$11]++;next}b[$11] < 2' - <(awk 'NR==FNR && !/^#/{a[$10]++;b[$11]++;next}a[$10] > 1 && a[$10] < 3' sample.intact.gff3 sample.intact.gff3) \
    | awk 'start == $4 || end == $5{print $0"\t"id}{start = $4; end = $5; id = $1"_"$4"_"$5"#"$3"\tlLTR";}' \
    | tr ' ' '\t' \
    | cut -f1-9,12,13 >> tmp_01

  ## rLTR
  awk 'NR==FNR && !/^#/{a[$10]++;b[$11]++;next}b[$11] > 1 && b[$11] < 3' sample.intact.gff3 sample.intact.gff3 \
    | awk 'NR==FNR && !/^#/{a[$10]++;b[$11]++;next}a[$10] < 2' - <(awk 'NR==FNR && !/^#/{a[$10]++;b[$11]++;next}b[$11] > 1 && b[$11] < 3' sample.intact.gff3 sample.intact.gff3) \
    | awk 'start == $4 || end == $5{print $0"\t"id}{start = $4; end = $5; id = $1"_"$4"_"$5"#"$3"\trLTR";}' \
    | tr ' ' '\t' \
    | cut -f1-9,12,13 >> tmp_01
  sort -k 4,4n -k 10,10 -k 11,11 tmp_01 > sample.intact.gff3

  ## remove complex structures; i.e., nested elements
  grep -v -f <(cut -f10 sample.intact.gff3 | var_count - | awk '!/^3/{split($2,arr,"#"); print arr[1]}') sample.intact.gff3 \
    | awk '$10 ~ $4 || $10 ~ $5' > tmp_01
  mv tmp_01 sample.intact.gff3

  ## number of raw_fl_LTRRT
  grep -c "fl_LTRRT" sample.intact.gff3 | tally_c T "raw_fl_LTRRT" - >> hit_count.txt
