## Align reads:
### Note: Cellranger v7+ counts intronic reads by default
### Note: using prebuilt reference

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/tran_2020/aligned_reads

cellranger="/opt/cellranger-7.1.0/bin/cellranger"

ref="/home/shared/hg_align_db/GRCh38_Cellranger_gencode_v32_ensembl_98/refdata-gex-GRCh38-2020-A"

for ea in $(ls ../raw_data); do
 $cellranger count \
  --id=$ea \
  --fastqs=../raw_data/$ea \
  --sample=$ea \
  --transcriptome=$ref \
  --localcores=10 \
  --no-bam
done

## Concatonate samples:

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/tran_2020/aligned_reads

cellranger="/opt/cellranger-7.1.0/bin/cellranger"

$cellranger aggr \
  --id="tran_2020_DLPFC" \
  --csv="cellranger_aggr.csv" \
  --normalize="none"
  
