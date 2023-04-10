## Get raw reads from SRA (control samples only):
## Acc. ID: SRP319543

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/morabito_2021

for ea in $(cat SRR_Acc_List.txt); do
  /opt/sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump -e 13 -m 100GB -S -p -O raw_data/ $ea
done

## Combine fastqs from the same sample

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/morabito_2021/raw_data

## GSM5292838: SRR14513977, SRR14513978, SRR14513979, SRR14513980, SRR14513981, SRR14513982, SRR14513983, SRR14513984
sampleid="GSM5292838"
filelist=$(ls | grep "$(seq 3977 3984)"[@])
mkdir $sampleid
cat $(echo "$filelist" | grep -E _2.fastq) > ${sampleid}/${sampleid}_S1_L001_R1_001.fastq
cat $(echo "$filelist" | grep -E _3.fastq) > ${sampleid}/${sampleid}_S1_L001_R2_001.fastq

## GSM5292851: SRR14514081, SRR14514082, SRR14514083, SRR14514084, SRR14514085, SRR14514086, SRR14514087, SRR14514088
sampleid="GSM5292851"
filelist=$(ls | grep "$(seq 4081 4088)"[@])
mkdir $sampleid
cat $(echo "$filelist" | grep -E _2.fastq) > ${sampleid}/${sampleid}_S1_L001_R1_001.fastq
cat $(echo "$filelist" | grep -E _3.fastq) > ${sampleid}/${sampleid}_S1_L001_R2_001.fastq

## GSM5292852: SRR14514089, SRR14514090, SRR14514091, SRR14514092, SRR14514093, SRR14514094, SRR14514095, SRR14514096
sampleid="GSM5292852"
filelist=$(ls | grep "$(seq 4089 4096)"[@])
mkdir $sampleid
cat $(echo "$filelist" | grep -E _2.fastq) > ${sampleid}/${sampleid}_S1_L001_R1_001.fastq
cat $(echo "$filelist" | grep -E _3.fastq) > ${sampleid}/${sampleid}_S1_L001_R2_001.fastq

## GSM5292853: SRR14514097, SRR14514098, SRR14514099, SRR14514100, SRR14514101, SRR14514102, SRR14514103, SRR14514104
sampleid="GSM5292853"
filelist=$(ls | grep "$(seq 4097 4104)"[@])
mkdir $sampleid
cat $(echo "$filelist" | grep -E _2.fastq) > ${sampleid}/${sampleid}_S1_L001_R1_001.fastq
cat $(echo "$filelist" | grep -E _3.fastq) > ${sampleid}/${sampleid}_S1_L001_R2_001.fastq

## GSM5292854: SRR14514105, SRR14514106, SRR14514107, SRR14514108, SRR14514109, SRR14514110, SRR14514111, SRR14514112
sampleid="GSM5292854"
filelist=$(ls | grep "$(seq 4105 4112)"[@])
mkdir $sampleid
cat $(echo "$filelist" | grep -E _2.fastq) > ${sampleid}/${sampleid}_S1_L001_R1_001.fastq
cat $(echo "$filelist" | grep -E _3.fastq) > ${sampleid}/${sampleid}_S1_L001_R2_001.fastq

## GSM5292855: SRR14514113, SRR14514114, SRR14514115, SRR14514116, SRR14514117, SRR14514118, SRR14514119, SRR14514120
sampleid="GSM5292855"
filelist=$(ls | grep "$(seq 4113 4120)"[@])
mkdir $sampleid
cat $(echo "$filelist" | grep -E _2.fastq) > ${sampleid}/${sampleid}_S1_L001_R1_001.fastq
cat $(echo "$filelist" | grep -E _3.fastq) > ${sampleid}/${sampleid}_S1_L001_R2_001.fastq

## GSM5292856: SRR14514121, SRR14514122, SRR14514123, SRR14514124, SRR14514125, SRR14514126, SRR14514127, SRR14514128
sampleid="GSM5292856"
filelist=$(ls | grep "$(seq 4121 4128)"[@])
mkdir $sampleid
cat $(echo "$filelist" | grep -E _2.fastq) > ${sampleid}/${sampleid}_S1_L001_R1_001.fastq
cat $(echo "$filelist" | grep -E _3.fastq) > ${sampleid}/${sampleid}_S1_L001_R2_001.fastq

## Align reads:
## Note: Cellranger v7+ counts intronic reads by default
## Note: using prebuilt reference

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/morabito_2021/aligned_reads

cellranger="/opt/cellranger-7.1.0/bin/cellranger"

ref="/home/shared/hg_align_db/GRCh38_Cellranger_gencode_v32_ensembl_98/refdata-gex-GRCh38-2020-A"

for ea in GSM5292854 GSM5292855 GSM5292856; do
  $cellranger count \
  --id=$ea \
  --fastqs=../raw_data/${ea} \
  --sample=$ea \
  --transcriptome=$ref \
  --localcores=10
  # --no-bam
done
