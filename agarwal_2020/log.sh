## Get raw reads from SRA (cortex samples only):
## Acc. ID: SRP229590

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/agarwal_2020/raw_data

for ea in $(cat ../SRR_Acc_List.txt); do 
  prefetch $ea --output-file ${ea}.sra
  fasterq-dump -b 100GB -c 100GB -m 100GB -S -p ${ea}.sra
done

## Edit filenames to be compatible with Cellranger:

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/agarwal_2020/raw_data

for ea in *_1.fastq; do mv $ea $(echo $ea | sed s/_1/_S1_L001_R1_001/); done
for ea in *_2.fastq; do mv $ea $(echo $ea | sed s/_2/_S1_L001_R2_001/); done

## Combine runs from the same sample:

n_files=$(cat ../rundata.csv | wc -w)

for i in $(seq 2 $n_files); do
  geo=$(awk -v i=$i -F',' 'NR==i{print $1}' ../rundata.csv)
  sra=$(awk -v i=$i -F',' 'NR==i{print $2}' ../rundata.csv)
  echo $geo
  if [[ ! -d $geo ]]; then mkdir $geo; fi
  ## Save file w/ list of all SRA runs associated with a sample (for grepping):
  awk -v i=$i -F',' 'NR==i{print $2; print $3}' ../rundata.csv | awk NF > runs.txt
  ## Grep all SRA runs:
  runsR1=($(printf '%s\n' * | grep -f runs.txt | grep '_R1_'))
  runsR2=($(printf '%s\n' * | grep -f runs.txt | grep '_R2_'))
  ## Concatonate runs from the same sample:
  cat ${runsR1[@]} > ${geo}/$(echo ${runsR1[1]} | sed s/$sra/$geo/)
  cat ${runsR2[@]} > ${geo}/$(echo ${runsR2[1]} | sed s/$sra/$geo/)
  ## Zip fastqs:
  pigz -p 5 ${geo}/*
done

# for ea in $samples; do mkdir $ea; mv ${ea}* $ea; done

## Align reads:
## Note: Cellranger v7+ counts intronic reads by default
## Note: using prebuilt reference

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/velmeshev_2019/aligned_reads

cellranger="/opt/cellranger-7.1.0/bin/cellranger"

ref="/home/shared/hg_align_db/GRCh38_Cellranger_gencode_v32_ensembl_98/refdata-gex-GRCh38-2020-A"

samples=$(awk -F',' 'NR>1{print substr($10,2,length($10)-2)}' ../sampleinfo.csv)

for ea in $samples; do
  $cellranger count \
    --id=$ea \
    --fastqs=../raw_data/${ea} \
    --sample=$ea \
    --transcriptome=$ref \
    --localcores=7
done