setwd("/mnt/bdata/rebecca/SCSN_meta_analysis/datasets/schirmer_2019")

library(dplyr)
library(data.table)

cellinfo <- fread("/mnt/bdata/@shared/scsn.expr_data/human_expr/postnatal/schirmer_2019/sampleinfo_barcode_level_tableS2.csv", data.table=F)

sampleinfo <- cellinfo %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(Batch=paste(unique(Capbatch), 
                               collapse=", "))

dat <- fread("SraRunTable.txt", data.table=F)
dat <-  dat %>%
  dplyr::select(Isolate, Run, Age, sex)

colnames(dat) <- c("Donor", "Accession_ID", "Age", "Sex")

dat <- merge(dat, sampleinfo, by=1)

write.table(dat, file="sampleinfo.csv", row.names=F, sep=",")
