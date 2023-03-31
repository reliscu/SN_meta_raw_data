## Get raw reads from SRA (control samples only):
### Acc. ID: SRP244344

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020

for ea in $(cat SRR_Acc_List.txt); do
  /home/shared/programs/sratoolkit.3.0.1-ubuntu64/bin/fasterq-dump -e 13 -m 100GB -S -p -O raw_data $ea
  for fastqs in $(raw_data/${ea}*); do
    pigz -p 13 raw_data/${fastqs}
  done
done

## Edit filenames to be compatible with Cellranger:
### E.g. sample_label_2.fastq.gz -> sample_label_S1_L001_R2_001.fastq.gz

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/raw_data

for ea in *_1.fastq; do
  mv $ea $(echo $ea | sed s/2/S1_L001_R1_001/) 
done

for ea in *_2.fastq; do
  mv $ea $(echo $ea | sed s/3/S1_L001_R2_001/)
done
   
## Align reads:
### Note: Cellranger v7+ counts intronic reads by default
### Note: using prebuilt reference

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/aligned_reads

cellranger="/opt/cellranger-7.1.0/bin/cellranger"

ref="/home/shared/hg_align_db/GRCh38_Cellranger_gencode_v32_ensembl_98/refdata-gex-GRCh38-2020-A"

for ea in $(cat ../SRR_Acc_List.txt); do
  $cellranger count \
    --id=$ea \
    --fastqs=../raw_data \
    --sample=$ea \
    --transcriptome=$ref \
    --localcores=14 \
    --no-bam
done

## Concatonate samples:

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/aligned_reads

cellranger="/opt/cellranger-7.1.0/bin/cellranger"

$cellranger aggr \
  --id="nagy_2020" \
  --csv="cellranger_aggr.csv" \
  --normalize="none"
  
