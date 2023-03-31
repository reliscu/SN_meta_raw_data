setwd("/mnt/bdata/rebecca/SCSN_meta_analysis/datasets/morabito_2021")

library(dplyr)
library(GEOquery)
library(data.table)

dat <- fread("SraRunTable.txt", data.table=F)
dat <- dat[grep("Control", dat$Diagnosis),]
dat <- dat[grep("nuclei", dat$isolation),]
dat <- dat[grep("TRANS", dat$LibrarySource),]

sampleinfo <- pData(getGEO("GSE174367")[[1]])
sampleinfo <- sampleinfo[is.element(sampleinfo$geo_accession, dat$`Sample Name`),]
sampleinfo[grep("101", sampleinfo$title), c(2)]
# [1] GSM5292839

## Sample 101 is not listed in paper supplementary material.
## According to GEO, it was not included in the analysis due to quality issues.

dat[is.element(dat$`Sample Name`, "GSM5292839"), c(1,2)]

## Manually removing runs associated with this sample from SRR accession list

dat <- dat[!is.element(dat$`Sample Name`, "GSM5292839"),]

## Get list of runs per sample:

dat %>% 
  dplyr::group_by(`Sample Name`) %>%
  dplyr::summarise(Runs=paste(unique(Run), collapse=", "))

