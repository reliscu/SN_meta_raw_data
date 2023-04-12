## Get raw reads from SRA:
## Acc. ID: SRP132816

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/velmeshev_2019/raw_data

for ea in $(cat ../SRR_Acc_List.txt); do
  nice -n -18 prefetch --output-directory $ea
done

for ea in SRR*sra; do
  run=$(echo $ea | sed 's/.sra//')
  nfastqs=$(printf '%s\n' ${run}* | grep fastq | wc -l)
  if [[ $nfastqs -lt 2 ]]; then
    echo $ea
    fasterq-dump -b 100GB -c 100GB -m 100GB -S -p $ea
  fi
done

## Edit filenames to be compatible with Cellranger:

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/velmeshev_2019/raw_data

for ea in *_1.fastq; do mv $ea $(echo $ea | sed s/_1/_S1_L001_R1_001/); done
for ea in *_2.fastq; do mv $ea $(echo $ea | sed s/_2/_S1_L001_R2_001/); done

## Make folders for each sample

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/velmeshev_2019/raw_data

samples=$(awk -F',' 'NR>1{print substr($10,2,length($10)-2)}' ../sampleinfo.csv)

for ea in $samples; do mkdir $ea; mv ${ea}* $ea; done

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