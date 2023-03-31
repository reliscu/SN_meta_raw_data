setwd("/mnt/bdata/rebecca/SCSN_meta_analysis/datasets/tran_2020/aligned_reads")

file_path <- "/mnt/bdata/rebecca/SCSN_meta_analysis/datasets/tran_2020/aligned_reads"

dat <- data.frame(sample_id=c("5161_DLPFC", "Br5207DLPFC", "BR5212_DLPFC"),
                  molecule_h5=c(file.path(file_path, "5161_DLPFC/outs/molecule_info.h5"),
                                file.path(file_path, "Br5207DLPFC/outs/molecule_info.h5"),
                                file.path(file_path, "BR5212_DLPFC/outs/molecule_info.h5")))

fwrite(dat, file="cellranger_aggr.csv")
