setwd("/mnt/bdata/rebecca/SCSN_meta_analysis/datasets/velmeshev_2019")

library(dplyr)
library(GEOquery)

sampleinfo <- read.csv("/home/shared/scsn.expr_data/human_expr/postnatal/velmeshev_2019/sampleinfo_sample_level_tableS1.csv")

dat <- fread("SraRunTable.txt", data.table=F)
dat <-  dat %>%
  dplyr::select(Run, `Library Name`)

sampleinfo <- sampleinfo %>%
  dplyr::filter(is.element(Sample.name, dat$`Library Name`)) %>%
  dplyr::select(Sample.name, Individual, Region, Age, 
                Sex, PMI, RIN, Diagnosis, Capbatch)

sampleinfo <- merge(sampleinfo, dat, by.x=1, by.y=2)
colnames(sampleinfo)[grep("Run", colnames(sampleinfo))] <- "Accession_ID"
colnames(sampleinfo)[grep("Capbatch", colnames(sampleinfo))] <- "Batch"
colnames(sampleinfo)[c(1)] <- "Label"

write.table(sampleinfo, file="sampleinfo.csv", row.names=F, sep=",")

## Add Cellranger cell estimates:

