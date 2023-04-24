## Align reads:
### Note: Cellranger v7+ counts intronic reads by default
### Note: using prebuilt reference

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/tran_2020/aligned_reads/temp

cellranger="/opt/cellranger-7.1.0/bin/cellranger"

ref="/home/shared/hg_align_db/GRCh38_Cellranger_gencode_v32_ensembl_98/refdata-gex-GRCh38-2020-A"

for ea in 5161_DLPFC  Br5207DLPFC  BR5212_DLPFC; do
 $cellranger count \
  --id=$ea \
  --fastqs=../../raw_data/$ea \
  --sample=$ea \
  --transcriptome=$ref \
  --localcores=10 
done
