setwd("/mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/aligned_reads")

file_path <- "/mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/aligned_reads"

sampleids <- unique(fread("../sampleinfo.csv", data.table=F)[,c(2)])

dat <- data.frame(sample_id=sampleids,
                  molecule_h5=file.path(file_path, sampleids, "outs/molecule_info.h5"))

fwrite(dat, file="cellranger_aggr.csv")
