## Get raw reads from SRA (control samples only):
## Acc. ID: SRP244344

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/raw_data

for ea in $(cat ../SRR_Acc_List.txt); do 
  prefetch $ea --output-file ${ea}.sra
  fasterq-dump -b 100GB -c 100GB -m 100GB -S -p ${ea}.sra
done

## Edit filenames to be compatible with Cellranger:

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/raw_data

for ea in *_1.fastq; do mv $ea $(echo $ea | sed s/2/S1_L001_R1_001/); done
for ea in *_2.fastq; do mv $ea $(echo $ea | sed s/3/S1_L001_R2_001/); done
   
## Combine and rename fastqs:

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/raw_data

nfiles=$(cat ../rundata.csv | wc -w)

for i in $(seq 2 $nfiles); do
  geo=$(awk -v i=$i -F',' 'NR==i{print $1}' ../rundata.csv)
  echo $geo
  if [[ ! -d $geo ]]; then mkdir $geo; fi
  ## Save file w/ list of all SRA runs associated with a sample (for grepping):
  awk -v i=$i -F',' 'NR==i{print $2; print $3}' ../rundata.csv | awk NF > runs.txt
  ## Grep all SRA runs:
  runsR1=($(printf '%s\n' * | grep -f runs.txt | grep '_R1_'))
  runsR2=($(printf '%s\n' * | grep -f runs.txt | grep '_R2_'))
  sra=$(awk -v i=$i -F',' 'NR==i{print $2}' ../rundata.csv)
  nruns=$(cat runs.txt | wc -w)
  if [[ $nruns -gt 1 ]]; then
    cat ${runsR1[@]} > ${geo}/$(echo ${runsR1[0]} | sed s/$sra/$geo/)
    cat ${runsR2[@]} > ${geo}/$(echo ${runsR2[0]} | sed s/$sra/$geo/)
  else
    ## Not all samples have multiple runs, and above code won't work in that case:
    cp $runsR1 ${geo}/$(echo $runsR1 | sed s/$sra/$geo/)
    cp $runsR2 ${geo}/$(echo $runsR2 | sed s/$sra/$geo/)
  fi
  ## Zip resultant fastqs:
  pigz -p 6 ${geo}/*
done

## Align reads:
## Note: Cellranger v7+ counts intronic reads by default
## Note: using prebuilt reference

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/aligned_reads

cellranger="/opt/cellranger-7.1.0/bin/cellranger"

ref="/home/shared/hg_align_db/GRCh38_Cellranger_gencode_v32_ensembl_98/refdata-gex-GRCh38-2020-A"

for ea in GSM*; do
  nice -n -18 $cellranger count \
    --id=$ea \
    --fastqs=../raw_data/${ea} \
    --sample=$ea \
    --transcriptome=$ref \
    --localcores=14
done
