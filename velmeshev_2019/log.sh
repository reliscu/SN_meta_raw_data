## Get raw reads from SRA:
## Acc. ID: SRP132816

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/velmeshev_2019/raw_data

for ea in $(cat ../SRR_Acc_List.txt); do
  nice -n -18 prefetch --output-directory $ea
done

for ea in SRR*; do fasterq-dump -m 100GB -S -p $ea; done
