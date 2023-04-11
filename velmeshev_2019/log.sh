## Get raw reads from SRA:
## Acc. ID: SRP132816

cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/velmeshev_2019/raw_data

for ea in $(cat ../SRR_Acc_List.txt); do
  nice -n -18 prefetch --output-directory $ea
done

for ea in SRR*sra; do
  run=$(echo $ea | sed 's/.sra//')
  nfastqs=$(printf '%s\n' ${run}* | grep fastq | wc -l)
  if [[ $nfastqs -lt 2 ]]; then
    echo $ea
    fasterq-dump -b 50GB -c 50GB -m 50GB -S -p $ea
  fi
done
