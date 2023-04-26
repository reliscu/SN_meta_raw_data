## Get raw reads from SRA (control samples only):
## Acc. ID: SRP199470 

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/schirmer_2019/raw_data

for ea in $(cat ../SRR_Acc_List.txt); do 
  prefetch $ea --max-size 200GB --output-file ${ea}.sra
  fasterq-dump -b 100GB -c 100GB -m 100GB -S -p ${ea}.sra
done

for ea in *fastq; do pigz -p 6 $ea; done

## Edit filenames to be compatible with Cellranger:

for ea in *_1.fastq.gz; do mv $ea $(echo $ea | sed s/_1/_S1_L001_R1_001/); done
for ea in *_2.fastq.gz; do mv $ea $(echo $ea | sed s/_2/_S1_L001_R2_001/); done

## Make folder per sample

for ea in $(cat ../SRR_Acc_List.txt); do 
  mkdir $ea
  mv ${ea}*.fastq.gz $ea
done

## Align reads:
## Note: Cellranger v7+ counts intronic reads by default
## Note: using prebuilt reference

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/schirmer_2019/aligned_reads

cellranger="/opt/cellranger-7.1.0/bin/cellranger"

ref="/home/shared/hg_align_db/GRCh38_Cellranger_gencode_v32_ensembl_98/refdata-gex-GRCh38-2020-A"

samples=$(awk -F',' 'NR>1{print substr($2,2,length($2)-2)}' ../sampleinfo.csv | uniq)

for ea in $samples; do
  $cellranger count \
    --id=$ea \
    --fastqs=../raw_data/${ea} \
    --sample=$ea \
    --transcriptome=$ref \
    --localcores=10
done
