setwd("/mnt/bdata/rebecca/SCSN_meta_analysis/datasets/nagy_2020/cluster_cellbender")

library(dplyr)
library(Seurat)

## Load data:

expr <- readRDS("../aligned_reads/nagy_2020/nagy_2020_cellbender.RDS")

dim(expr)
# [1] 36601 29809

## QC already performed:

VlnPlot(expr, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

summary(expr[[]]$nCount_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 547    1145    2394    3786    4893   49403 
summary(expr[[]]$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 501     868    1604    1988    2725    8968  

expr[["Sample"]] <- as.character(sapply(strsplit(colnames(expr), "-", fixed=T), "[", 2))

## Cells per sample:

table(expr$Sample)
# 1   10   11   12   17   18   19    2    3    4    5    6    7    8 
# 2019 4103 1930 2175 1061 2292 3190 1179 2344 1205 1714 1728 3589 1280 

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
expr <- IntegrateData(anchorset=anchors, normalization.method="SCT")
expr <- RunPCA(expr, npcs=20)
expr <- RunUMAP(expr, reduction="pca", dims=1:20, min.dist=.5)
expr <- FindNeighbors(expr, reduction="pca", dims=1:20)
expr <- FindClusters(expr, resolution=.1)

pdf("initial_clusters_cellbender.pdf")
DimPlot(expr, reduction="umap", pt.size=.3, label=T)
DimPlot(expr, reduction="umap", group.by="Sample", pt.size=.3, label=T)
FeaturePlot(expr, features="nFeature_RNA") 
FeaturePlot(expr, features="percent.mt")
dev.off()

## Visualize canonical markers:

DefaultAssay(object=expr) <- "SCT"

pdf("canonical_markers_cellbender.pdf")

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
expr$Cell_Class[is.element(expr$Cell_Class, c(6))] <- "OG"
expr$Cell_Class[is.element(expr$Cell_Class, c(8))] <- "OPC"
expr$Cell_Class[is.element(expr$Cell_Class, c(1))] <- "ASC"
expr$Cell_Class[is.element(expr$Cell_Class, c(0, 2, 3, 7, 11:13, 15, 16))] <- "EXC"
expr$Cell_Class[is.element(expr$Cell_Class, c(4, 5, 9, 10, 14))] <- "INH"

Idents(expr) <- expr$Cell_Class

pdf("annotated_clusters_cellbender.pdf")
DimPlot(expr, reduction="umap", pt.size=.3, label=T)
dev.off()

saveRDS(expr, file="expression_nagy_2020_cellbender_annotated.RDS")