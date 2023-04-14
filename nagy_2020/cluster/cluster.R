setwd("/mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/cluster")

library(dplyr)
library(Seurat)

## Load data:

expr <- readRDS("../aligned_reads/nagy_2020/nagy_2020.RDS")
dim(expr)
# [1] 36601 33260

## QC already performed:

VlnPlot(expr, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

summary(expr[[]]$nCount_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 552    1180    2368    3716    4739   49543 
summary(expr[[]]$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 501     900    1599    1980    2680    8995 

expr[["Sample"]] <- as.character(sapply(strsplit(colnames(expr), "_", fixed=T), "[", 1))

## Cells per sample:

table(expr$Sample)
# GSM4281975 GSM4281976 GSM4281980 GSM4281982 GSM4281985 GSM4281986 GSM4281988 
# 2013       1153       2340       1204       1691       1727       3541 
# GSM4281989 GSM4281992 GSM4281993 GSM4281994 GSM4281995 GSM4281997 GSM4281998 
# 1253       2750       4037       1920       2142        940         76 
# GSM4282000 GSM4282002 GSM4282004 
# 1067       2291       3115 

## Integrate samples:
### Ref: https://satijalab.org/seurat/articles/integration_introduction.html

split <- SplitObject(expr, split.by="Sample")
split <- lapply(X=split, FUN=function(x){
  set.seed(123)
  x <- SCTransform(x, vars.to.regress=c('percent.mt', 'nCount_RNA'), vst.flavor="v2")
})
features <- SelectIntegrationFeatures(object.list=split, nfeatures=2000)
split <- PrepSCTIntegration(object.list=split, anchor.features=features)
anchors <- FindIntegrationAnchors(object.list=split, anchor.features=features, normalization.method="SCT")
kw <- 100
if(min(table(expr$Sample))<100){
  kw <- min(table(expr$Sample))
}
expr <- IntegrateData(anchorset=anchors, normalization.method="SCT", k.weight=kw)
expr <- RunPCA(expr, npcs=20)
expr <- RunUMAP(expr, reduction="pca", dims=1:20, min.dist=.5)
expr <- FindNeighbors(expr, reduction="pca", dims=1:20)
expr <- FindClusters(expr, resolution=.1)

pdf("initial_clusters.pdf")
DimPlot(expr, reduction="umap", pt.size=.3, label=T)
DimPlot(expr, reduction="umap", group.by="Sample", pt.size=.3, label=T)
FeaturePlot(expr, features="nFeature_RNA") 
FeaturePlot(expr, features="percent.mt")
dev.off()

## Visualize canonical markers:

DefaultAssay(object=expr) <- "SCT"

pdf("canonical_markers.pdf")

## OG
genes <- c("MAG", "MOG")
FeaturePlot(object=expr, features=genes, label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)

## OPC
genes <- c("PDGFRA", "PCDH15", "VCAN")
FeaturePlot(object=expr, features=genes, label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)

## ASC
genes <- c("GFAP", "SLC1A2", "AQP4")
FeaturePlot(object=expr, features=genes, label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)

## EXC
genes <- c("SLC17A7", "SLC1A1", "SLC17A6")
FeaturePlot(object=expr, features=genes, label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)

## INH
genes <- c("GAD2", "GAD2", "SLC6A1")
FeaturePlot(object=expr, features=genes, label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)

# END
genes <- c("PECAM1", "ID3", "FLT1", "TNFRSF21")
FeaturePlot(object=expr, features=genes, label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)

# PER
genes <- c("ACTA2", "PDGFRB", "MCAM")
FeaturePlot(object=expr, features=genes, label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)

# EPEN
genes <- c("FOXJ1")
FeaturePlot(object=expr, features=genes, label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)

# FIB
genes <- c("COL1A1", "LUM", "COL6A2")
FeaturePlot(object=expr, features=genes, label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)

## MIC
genes <- c("ITGAM", "AIF1", "TMEM119")
FeaturePlot(object=expr, features=genes, label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)

dev.off()

## Find cluster markers:

expr <- PrepSCTFindMarkers(expr)
markers <- FindAllMarkers(expr, assay="SCT")

markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::slice_head(n=10) %>%
  as.data.frame()

expr$Cell_Class <- as.character(expr$seurat_clusters)

expr$Cell_Class[is.element(expr$Cell_Class, c(0, 1))] <- "OG"
expr$Cell_Class[is.element(expr$Cell_Class, c(5))] <- "OPC"
expr$Cell_Class[is.element(expr$Cell_Class, c(2))] <- "ASC"
expr$Cell_Class[is.element(expr$Cell_Class, c(4))] <- "EXC"
expr$Cell_Class[is.element(expr$Cell_Class, c(6, 7))] <- "INH"
expr$Cell_Class[is.element(expr$Cell_Class, c(3))] <- "MIC"
expr$Cell_Class[is.element(expr$Cell_Class, c(9))] <- "PER"
expr$Cell_Class[is.element(expr$Cell_Class, c(8, 10))] <- "NEU"

Idents(expr) <- expr$Cell_Class

pdf("annotated_clusters.pdf")
DimPlot(expr, reduction="umap", pt.size=.3, label=T)
dev.off()

DefaultAssay(object=expr) <- "RNA"

saveRDS(expr, file="expression_nagy_2020_annotated.RDS")
