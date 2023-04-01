setwd("/mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020")

library(dplyr)
library(GEOquery)

sampleinfo <- pData(getGEO("GSE174367")[[1]])

dat <- fread("SraRunTable.txt", data.table=F)

sampleinfo <- sampleinfo %>%
  dplyr::filter(is.element(geo_accession, sampleids)) %>%
  dplyr::select(title, 
                geo_accession, 
                `age:ch1`, 
                `Sex:ch1`)

colnames(sampleinfo) <- c("Label", "Accession_ID", "Age", "Sex")

## Per table S1:

sampleinfo$Batch <- 2
sampleinfo$Batch[grep("100|96", sampleinfo$Label)] <- 3

## Add Cellranger cell estimates:

sampleinfo$No.Nuclei <- 4600
sampleinfo$No.Nuclei[sampleinfo$Accession_ID=="GSM5292851"] <- 2800
sampleinfo$No.Nuclei[sampleinfo$Accession_ID=="GSM5292852"] <- 2600
sampleinfo$No.Nuclei[sampleinfo$Accession_ID=="GSM5292853"] <- 2800
sampleinfo$No.Nuclei[sampleinfo$Accession_ID=="GSM5292854"] <- 2400
sampleinfo$No.Nuclei[sampleinfo$Accession_ID=="GSM5292855"] <- 2400
sampleinfo$No.Nuclei[sampleinfo$Accession_ID=="GSM5292856"] <- 4500

write.table(sampleinfo, file="sampleinfo.csv", row.names=F, sep=",")
