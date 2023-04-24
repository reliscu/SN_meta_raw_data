setwd("/mnt/bdata/rebecca/SCSN_meta_analysis/datasets/agarwal_2020")

library(dplyr)
library(GEOquery)
library(data.table)

sampleinfo <- pData(getGEO("GSE140231")[[1]])

sampleinfo <- sampleinfo %>%
  dplyr::select(`title`, `geo_accession`, `age:ch1`, `batch:ch1`, 
                `gender:ch1`, `post-mortem-interval(hours):ch1`,)

colnames(sampleinfo) <- c("Label", "Accession_ID", "Age", "Batch", "Sex", "PMI")

dat <- fread("SraRunTable.txt", data.table=F)
dat <-  dat %>%
  dplyr::select(Run, `Sample Name`)

sampleinfo <- merge(sampleinfo, dat, by=2)

rundata <- sampleinfo %>%
  dplyr::group_by(Accession_ID) %>%
  dplyr::summarise(Run=paste(unique(Run), collapse=","))

write.table(sampleinfo, file="sampleinfo.csv", row.names=F, sep=",")
write.table(rundata, file="rundata.csv", row.names=F, sep=",", quote=F)
