#!/usr/bin/env bash

## Combine and rename fastqs

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/raw_data

n_files=$(cat ../rundata.csv | wc -w)

for i in $(seq 2 $n_files); do
  geo=$(awk -v i=$i -F',' 'NR==i{print $1}' ../rundata.csv)
  sra=$(awk -v i=$i -F',' 'NR==i{print $2}' ../rundata.csv)
  
  if [[ ! -d $geo ]]; then mkdir $geo; fi

  ## Save file with all SRA runs associated with a sample:
  awk -v i=$i -F',' 'NR==i{print $2; print $3}' ../rundata.csv | awk NF > runs.txt
  runsR1=($(printf '%s\n' * | grep -f runs.txt | grep '_R1_'))
  runsR2=($(printf '%s\n' * | grep -f runs.txt | grep '_R2_'))
  
  n_runs=$(cat runs.txt | wc -w)
  if [[ $n_runs > 1 ]] 
  then
    cat ${runsR1[@]} > ${geo}/$(echo ${runsR2[1]} | sed s/$sra/$geo/)
    cat ${runsR1[@]} > ${geo}/$(echo ${runsR2[1]} | sed s/$sra/$geo/)
  else
    cp ${runsR1[@]} ${geo}/$(echo $runsR1 | sed s/$sra/$geo/)
    cp ${runsR2[@]} ${geo}/$(echo $runsR2 | sed s/$sra/$geo/)
  fi
done

rm runs.txt

## Align reads:
### Note: Cellranger v7+ counts intronic reads by default
### Note: using prebuilt reference

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/aligned_reads

cellranger="/opt/cellranger-7.1.0/bin/cellranger"

ref="/home/shared/hg_align_db/GRCh38_Cellranger_gencode_v32_ensembl_98/refdata-gex-GRCh38-2020-A"

samples=$(awk -F',' 'NR>1{print $1}' ../rundata.csv)

for ea in $samples; do
  $cellranger count \
  --id=$ea \
  --fastqs=../raw_data/${ea} \
  --sample=$ea \
  --transcriptome=$ref \
  --localcores=14
done