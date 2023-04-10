setwd("/mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020")

library(dplyr)
library(GEOquery)

sampleinfo <- pData(getGEO("GSE144136")[[1]])

dat <- read.csv("/mnt/bdata/@shared/scsn.expr_data/human_expr/postnatal/nagy_2020/sampleinfo_batch_from_authors.csv")
dat <- dat[grep("Control", dat$Condition),]
table(dat$Batch)
# B1 B2 B3 B4 B5 B6 
# 1  1  4  4  4  5 

sampleinfo <- sampleinfo %>%
  dplyr::select(title, 
                geo_accession, 
                `Sex:ch1`,
                `group:ch1`) %>%
  dplyr::mutate(
    title=as.character(
      sapply(strsplit(title, ":"), "[", 1)
    )
  )

colnames(sampleinfo) <- c("Label", "Accession_ID", "Sex", "Diagnosis")

sampleinfo <- merge(sampleinfo, dat, by=1)
sampleinfo %>%
  dplyr::group_by(Batch) %>%
  dplyr::summarise(n=paste(sort(Sex), collapse=", "))

runtable <- read.csv("SraRunTable.txt")
runtable <- runtable %>%
  dplyr::select(Run, Sample.Name)

sampleinfo <- merge(sampleinfo, runtable, by.x="Accession_ID", by.y="Sample.Name")

rundata <- sampleinfo %>%
  dplyr::group_by(Accession_ID) %>%
  dplyr::summarise(Run=paste(unique(Run), collapse=","))

write.table(sampleinfo, file="sampleinfo.csv", row.names=F, sep=",")
write.table(rundata, file="rundata.csv", row.names=F, sep=",", quote=F)

# sampleinfo <- read.csv("sampleinfo.csv")


