## Get raw reads from SRA (control samples only):
## Acc. ID: SRP244344

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020

for ea in $(cat SRR_Acc_List.txt); do
  /home/shared/programs/sratoolkit.3.0.1-ubuntu64/bin/fasterq-dump -e 13 -m 100GB -S -p -O raw_data $ea
  for fastqs in $(raw_data/${ea}*); do
    pigz -p 13 raw_data/${fastqs}
  done
done

## Edit filenames to be compatible with Cellranger:

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/raw_data

for ea in *_1.fastq; do mv $ea $(echo $ea | sed s/2/S1_L001_R1_001/); done
for ea in *_2.fastq; do mv $ea $(echo $ea | sed s/3/S1_L001_R2_001/); done
   
## Combine and rename fastqs:

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/raw_data

n_files=$(cat ../rundata.csv | wc -w)

for i in $(seq 2 $n_files); do
  geo=$(awk -v i=$i -F',' 'NR==i{print $1}' ../rundata.csv)
  sra=$(awk -v i=$i -F',' 'NR==i{print $2}' ../rundata.csv)
  
  if [[ ! -d $geo ]]; then mkdir $geo; fi

  echo $geo
  
  ## Save file w/ list of all SRA runs associated with a sample (for grepping):
  awk -v i=$i -F',' 'NR==i{print $2; print $3}' ../rundata.csv | awk NF > runs.txt
  ## Grep all SRA runs:
  runsR1=($(printf '%s\n' * | grep -f runs.txt | grep '_R1_'))
  runsR2=($(printf '%s\n' * | grep -f runs.txt | grep '_R2_'))
  
  n_runs=$(cat runs.txt | wc -w)
  if [[ $n_runs > 1 ]]; then
  
    # cat ${runsR1[@]} > ${geo}/$(echo ${runsR1[1]} | sed s/$sra/$geo/)
    # cat ${runsR2[@]} > ${geo}/$(echo ${runsR2[1]} | sed s/$sra/$geo/)
    cp ${runsR1[@]} ${geo}
    cp ${runsR2[@]} ${geo}
    
  else
    cp $runsR1 ${geo}/$(echo $runsR1 | sed s/$sra/$geo/)
    cp $runsR2 ${geo}/$(echo $runsR2 | sed s/$sra/$geo/)
  fi
done

rm runs.txt

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/raw_data/GSM4281997
for ea in SRR12391931*; do mv $ea $(echo $ea | sed s/S1/S2/); done
for ea in *; do mv $ea $(echo $ea | sed -E s/'SRR[0-9]+'/GSM4281997/); done

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/raw_data/GSM4281998
for ea in SRR12391933*; do mv $ea $(echo $ea | sed s/S1/S2/); done
for ea in *; do mv $ea $(echo $ea | sed -E s/'SRR[0-9]+'/GSM4281998/); done

## Align reads:
## Note: Cellranger v7+ counts intronic reads by default
## Note: using prebuilt reference

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/aligned_reads

cellranger="/opt/cellranger-7.1.0/bin/cellranger"

ref="/home/shared/hg_align_db/GRCh38_Cellranger_gencode_v32_ensembl_98/refdata-gex-GRCh38-2020-A"

samples=$(awk -F',' 'NR>1{print $1}' ../rundata.csv)

for ea in GSM4281997 GSM4281998; do
  nice -n -18 $cellranger count \
    --id=$ea \
    --fastqs=../raw_data/${ea} \
    --sample=$ea \
    --transcriptome=$ref \
    --localcores=14
done

## Concatonate samples:

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/aligned_reads

cellranger="/opt/cellranger-7.1.0/bin/cellranger"

$cellranger aggr \
  --id="nagy_2020" \
  --csv="cellranger_aggr.csv" \
  --normalize="none"
  
