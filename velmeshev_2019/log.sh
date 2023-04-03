## Get raw reads from SRA (control samples only):
## Acc. ID: SRP132816

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/velmeshev_2019

for ea in $(cat SRR_Acc_List.txt); do
  /home/shared/programs/sratoolkit.3.0.1-ubuntu64/bin/fasterq-dump -e 13 -m 100GB -S -p -O raw_data $ea
  for fastqs in $(raw_data/${ea}*); do
    pigz -p 13 raw_data/${fastqs}
  done
done