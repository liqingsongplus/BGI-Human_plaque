## Combined with the single cell annotation results, the spatial group data were annotated.

## The Stereo-seq data could be obtained from the CNGB Spatial Transcript Omics 
## DataBase (https://db.cngb.org/stomics/) with accession number: STT0000052.

library(Seurat)
library(tidyverse)

setwd("D:/work/3_4_human_plaque_2024/1_raw_sp")  # Set the Stereo-seq data path

sc = readRDS("D:/work/3_4_human_plaque_2024/11_norm_PLA/HP_harmony_merge.rds") # Read integrated single cell data.

DefaultAssay(sc) = 'RNA'
sc = sc |> 
  SCTransform(assay = "RNA") %>% # , verbose = FALSE
  RunPCA() %>%
  RunUMAP(dims = 1:30)


A = dir()  # Get all Stereo-seq names
for(file in A){
  print(file)
  seu = readRDS(file)
  seu <- SCTransform(seu, assay = "Spatial") %>% 
    RunPCA()%>%
    RunUMAP(dims = 1:30)
  
  anchors <- FindTransferAnchors(reference = sc, query = seu, normalization.method = "SCT")
  
  predictions.assay <- TransferData(anchorset = anchors, refdata = sc$sub_cell, 
                                    weight.reduction = seu[["pca"]], dims = 1:30) 
  
  pre = gsub("_filter.rds", "", file)
  saveRDS(predictions.assay, paste0(pre,"_predicted_assay.rds"))
  
  celltype_predict = predictions.assay[c('predicted.id','prediction.score.max')]
  seu <- AddMetaData(seu, metadata = celltype_predict)
  table(seu$predicted.id)
  
  saveRDS(seu@meta.data,  paste0(pre,"_meta.rds"))  # Save stereo-seq annotation data into metadata.
}
