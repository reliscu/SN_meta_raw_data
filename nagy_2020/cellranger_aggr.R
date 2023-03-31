setwd("/mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/aligned_reads")

file_path <- "/mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/aligned_reads"

sample_ids <- fread("../SRR_Acc_List.txt", header=F, data.table=F)[,c(1)]

dat <- data.frame(sample_id=sample_ids,
                  molecule_h5=file.path(file_path, sample_ids, "outs/molecule_info.h5"))

fwrite(dat, file="cellranger_aggr.csv")
