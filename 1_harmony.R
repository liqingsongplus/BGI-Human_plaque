## After single cell quality control, integration, clustering and batch removal were carried out
## The integrated single-cell data can be downloaded from CNSA using accession number: CNP0004894. (https://db.cngb.org/cnsa)


library(harmony)
library(Seurat)
library(tidyverse)

setwd("D:/work/3_4_human_plaque_2024/11_norm_PLA/")

seu = NULL

A = dir()
for(i in A){
  if(is.null(seu)){
    seu = readRDS(paste0(i,"/",i, "_QC.rds")) # 
  }else{
    y = readRDS(paste0(i,"/",i, "_QC.rds"))
    y = subset(y, nFeature_RNA > 500)
    seu = merge(seu, y)
    print(i)
    rm(y)
    gc()
  }
}

seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
seu <- harmony::RunHarmony(seu, group.by.vars = "orig.ident")  # Remove batch effect by sample source
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30 )
seu = FindClusters(seu,resolution = c(0.2,0.4,0.6,0.8,1,1.2)) 
saveRDS(seu,"HP_harmony_merge.rds")


