cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/aligned_reads

source activate CellBender

for ea in GSM*; do
  if [[ ! -d ${ea}/cellbender ]]; then
    mkdir ${ea}/cellbender
  fi
  ncells=$(awk -F',' 'NR>1{print substr($1$2,2,length($1$2)-2)}' ${ea}/outs/metrics_summary.csv)
  nice -n 10 cellbender remove-background \
    --input ${ea}/outs/raw_feature_bc_matrix.h5 \
    --output ${ea}/cellbender/output.h5 \
    --cuda \
    --expected-cells $ncells \
    --total-droplets-included 10000
done
