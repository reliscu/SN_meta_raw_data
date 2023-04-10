## Get raw reads from SRA (control samples only):
## Acc. ID: SRP132816

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/velmeshev_2019

fasterq="/home/shared/programs/sratoolkit.3.0.1-ubuntu64/bin/fasterq-dump"

for ea in $(cat SRR_Acc_List.txt); do
  nice -n -18 prefetch --output-directory raw_data $ea
 # nice -n -18 $fasterq -m 100GB -S -p -O raw_data/${ea}
done
