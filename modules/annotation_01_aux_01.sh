#!/bin/bash

set -e
set -u
set -o pipefail

  ## Upstream/downstream sequences of k bp length are locally aligned to locate TSDs
  ## by default, 6 bp are surveyed and up to 1 mismatch is allowed
  ## parrallel run

  ## load lib
  source ${bin_path}/pfunlib.sh

  anno_batch_=${1}

  tmpdir_=$(mktemp -d tmp.XXXXXXXXXX)
  outdir_=$PWD

  ## retrieve TSD seq case by case, and annotate
  test_rm_ith "$tmpdir_/soloLTR.gff3"
  counter=0
  awk '{split($4,arr,"#"); print arr[1]}' $anno_batch_ | while read -r solo; do

    counter=$[$counter +1]

    counter_id=$(wc -l < $anno_batch_ | awk '{print length($0)}' | awk -v counter=${counter} '{pad = $0; printf "%0*d", pad, counter}')
    strand=$(awk -v solo=${solo} '$4 ~ "_5p" && $4 ~ "#"solo"_"{print $6}' LTR_blast_single_hit_int_dom.bed)

    mismatch_n=$(awk -v solo=${solo} '$2 ~ "#"solo{print $8}' ltr_tsd)
    alg_length=$(awk -v solo=${solo} '$2 ~ "#"solo{print $7}' ltr_tsd)

    motif=$(tr '\n' ' ' < sample.intact.fa--soloLTR | sed 's/>/\n/g' | awk -v solo=${solo} '$1 ~ solo{print substr($2,1,2)substr($2,length($2)-1,2)}')

    ## TSD seq for soloLTR(-) are not RC
    if [ "${strand}" == "minus" ]; then
      ## following condition means that the mismatch must be located at either termini; that's evaluated in $tsd_v2_test
      if [ "${alg_length}" -eq 4 ] && [ "${mismatch_n}" -eq 0 ]; then

        ltsd=$(grep -A1 "${solo}.*5p" LTR_blast_single_hit.fa \
                 | tail -n1 \
                 | awk '{print toupper(substr($0,2,5))}' \
                 | reverse_complement -)
        rtsd=$(grep -A1 "${solo}.*3p" LTR_blast_single_hit.fa \
                 | tail -n1 \
                 | awk '{print toupper(substr($0,1,5))}' \
                 | reverse_complement -)
        tsd_v2_test=$(paste -d '\t' <(echo ${ltsd} | grep -o '.') <(echo ${rtsd} | grep -o '.') \
                        | awk '$1==$2' \
                        | wc -l \
                        | awk '{if($0<4){print "ERR"}else{print "TP"}}')

      else

        ltsd=$(grep -P "#${solo}\t" ltr_tsd \
                 | cut -f3 \
                 | awk '{print toupper($0)}' \
                 | reverse_complement -)
        rtsd=$(grep -P "#${solo}\t" ltr_tsd \
                 | cut -f4 \
                 | awk '{print toupper($0)}' \
                 | reverse_complement -)
        tsd_v2_test="TP"

        if [ "${alg_length}" -eq 6 ] && [ "${mismatch_n}" -gt 1 ]; then
          tsd_v2_test="ERR"
        fi

      fi
    else
      ## following condition means that the mismatch must be located at either termini; that's evaluated in $tsd_v2_test
      if [ "${alg_length}" -eq 4 ] && [ "${mismatch_n}" -eq 0 ]; then

        ltsd=$(grep -A1 "${solo}.*5p" LTR_blast_single_hit.fa | tail -n1 | awk '{print substr($0,2,5)}')
        rtsd=$(grep -A1 "${solo}.*3p" LTR_blast_single_hit.fa | tail -n1 | awk '{print substr($0,1,5)}')

        tsd_v2_test=$(paste -d '\t' <(echo ${ltsd} | grep -o '.') <(echo ${rtsd} | grep -o '.') \
                        | awk '$1==$2' \
                        | wc -l \
                        | awk '{if($0<4){print "ERR"}else{print "TP"}}')

      else

        ltsd=$(grep -P "#${solo}\t" ltr_tsd | cut -f3)
        rtsd=$(grep -P "#${solo}\t" ltr_tsd | cut -f4)
        tsd_v2_test="TP"

        if [ "${alg_length}" -eq 6 ] && [ "${mismatch_n}" -gt 1 ]; then
          tsd_v2_test="ERR"
        fi
      fi
    fi


    ## generate gff3
    if [ "${tsd_v2_test}" != "ERR" ]; then
      awk -v solo=${solo} '$4 ~ solo"#"' $anno_batch_ \
        | awk -v sample=${sample} -v ltsd=${ltsd} -v rtsd=${rtsd} -v motif=${motif} -v mismatch_n=${mismatch_n} -v strand=${strand} -v counter_id=${counter_id} -v alg_length=${alg_length} 'BEGIN{OFS="\t"}{if(mismatch_n == 0 && alg_length == 4){split($4,arr,"#|/"); print $1, "soloLTRseeker", "soloLTR", $2+1, $3, ".", strand, ".", "ID=soloLTR/"sample"/"$1"/"$2+1"/"$3";Name=soloLTR_"counter_id";Classification=LTR/"arr[2]";Sequence_ontology=SO:0001003;Lineage="arr[3]";motif="toupper(motif)";TSD="toupper(ltsd)"_"toupper(rtsd)"_1"}else{split($4,arr,"#|/"); print $1, "soloLTRseeker", "soloLTR", $2+1, $3, ".", strand, ".", "ID=soloLTR/"sample"/"$1"/"$2+1"/"$3";Name=soloLTR_"counter_id";Classification=LTR/"arr[2]";Sequence_ontology=SO:0001003;Lineage="arr[3]";motif="toupper(motif)";TSD="toupper(ltsd)"_"toupper(rtsd)"_"mismatch_n}}' >> $tmpdir_/soloLTR.gff3
    fi

  done
