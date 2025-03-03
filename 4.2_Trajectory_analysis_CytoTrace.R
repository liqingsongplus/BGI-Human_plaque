## A pseudo-time sequence analysis of a given cell type was performed using CytoTRACE

library(CytoTRACE)
library(Seurat)
library(tidyverse)

scRNA <- readRDS('B_sub.rds')  # Import the single-cell RNA sequencing data stored in Seurat object format to perform trajectory analysis.
result <- CytoTRACE(as.matrix(scRNA@assays$RNA@counts),ncores = 1)
saveRDS(result,file = 'B_sub_cyto.rds')

pos <- names(head(result$cytoGenes,1500)) #top postive 1500
neg <- names(tail(result$cytoGenes,1500)) #top negtive 1500
scRNA <- NormalizeData(scRNA,assay = 'RNA')
var_gene <- c(pos,neg)
VariableFeatures(scRNA) <- var_gene

DefaultAssay(scRNA) = "RNA"
phe <- scRNA$sub_cell  # Group according to the specified cell type
phe <- as.character(phe)
names(phe) <- rownames(scRNA@meta.data)
plotCytoTRACE(result,phenotype = phe,emb = scRNA@reductions$umap@cell.embeddings)
plotCytoGenes(result, numOfGenes = 10)
write.csv(result$cytoGenes,'CytoTRACE_genes.csv')
