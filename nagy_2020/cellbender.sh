cd /mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/aligned_reads

source activate CellBender

mkdir SRR12391908/cellbender
cellbender remove-background \
  --input SRR12391908/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391908/output.h5 \
  --cuda \
  --expected-cells 2000 \
  --total-droplets-included 20000

mkdir SRR12391909/cellbender
cellbender remove-background \
  --input SRR12391909/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391909/cellbender/output.h5 \
  --cuda \
  --expected-cells 1500 \
  --total-droplets-included 20000

mkdir SRR12391913/cellbender
cellbender remove-background \
  --input SRR12391913/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391913/cellbender/output.h5 \
  --cuda \
  --expected-cells 2000 \
  --total-droplets-included 20000

mkdir SRR12391915/cellbender
cellbender remove-background \
  --input SRR12391915/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391915/cellbender/output.h5 \
  --cuda \
  --expected-cells 1200 \
  --total-droplets-included 20000

mkdir SRR12391918/cellbender
cellbender remove-background \
  --input SRR12391918/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391918/cellbender/output.h5 \
  --cuda \
  --expected-cells 2600 \
  --total-droplets-included 20000

mkdir SRR12391919/cellbender
cellbender remove-background \
  --input SRR12391919/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391919/cellbender/output.h5 \
  --cuda \
  --expected-cells 2300 \
  --total-droplets-included 20000

mkdir SRR12391921/cellbender
cellbender remove-background \
  --input SRR12391921/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391921/cellbender/output.h5 \
  --cuda \
  --expected-cells 3800 \
  --total-droplets-included 20000

mkdir SRR12391922/cellbender
cellbender remove-background \
  --input SRR12391922/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391922/cellbender/output.h5 \
  --cuda \
  --expected-cells 1200 \
  --total-droplets-included 20000

#### START HERE #####

mkdir SRR12391925/cellbender
cellbender remove-background \
  --input SRR12391925/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391925/cellbeder/output.h5 \
  --cuda \
  --expected-cells 2500 \
  --total-droplets-included 20000

mkdir SRR12391926/cellbender
cellbender remove-background \
  --input SRR12391926/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391926/cellbender/output.h5 \
  --cuda \
  --expected-cells 4400 \
  --total-droplets-included 20000

mkdir SRR12391927/cellbender
cellbender remove-background \
  --input SRR12391927/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391927/cellbender/output.h5 \
  --cuda \
  --expected-cells 2200 \
  --total-droplets-included 20000
 
mkdir SRR12391928/cellbender 
cellbender remove-background \
  --input SRR12391928/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391928/cellbender/output.h5 \
  --cuda \
  --expected-cells 2000 \
  --total-droplets-included 20000

mkdir SRR12391935/cellbender   
cellbender remove-background \
  --input SRR12391935/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391935/cellbender/output.h5 \
  --cuda \
  --expected-cells 1000 \
  --total-droplets-included 20000

mkdir SRR12391937/cellbender    
cellbender remove-background \
  --input SRR12391937/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391937/cellbender/output.h5 \
  --cuda \
  --expected-cells 2000 \
  --total-droplets-included 20000

mkdir SRR12391939/cellbender
cellbender remove-background \
  --input SRR12391939/outs/raw_feature_bc_matrix.h5 \
  --output SRR12391939/cellbender/output.h5 \
  --cuda \
  --expected-cells 3000 \
  --total-droplets-included 20000