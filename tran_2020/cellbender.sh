cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/tran_2020/aligned_reads

source activate CellBender

mkdir 5161_DLPFC/cellbender
cellbender remove-background \
  --input 5161_DLPFC/outs/raw_feature_bc_matrix.h5 \
  --output 5161_DLPFC/cellbender/output.h5 \
  --cuda \
  --expected-cells 4500 \
  --total-droplets-included 10000

mkdir Br5207DLPFC/cellbender
cellbender remove-background \
  --input Br5207DLPFC/outs/raw_feature_bc_matrix.h5 \
  --output Br5207DLPFC/cellbender/output.h5 \
  --cuda \
  --expected-cells 6000 \
  --total-droplets-included 10000

mkdir BR5212_DLPFC/cellbender
cellbender remove-background \
  --input BR5212_DLPFC/outs/raw_feature_bc_matrix.h5 \
  --output BR5212_DLPFC/cellbender/output.h5 \
  --cuda \
  --expected-cells 1800 \
  --total-droplets-included 10000
