## A pseudo-time sequence analysis of a given cell type was performed using monocle3

library(monocle3)
library(Matrix)
library(Seurat)
library(ggplot2)

seu = readRDS("B_sub.rds")  # Import the single-cell RNA sequencing data stored in Seurat object format to perform trajectory analysis.
Idents(seu) =seu$sub_cell

data <- GetAssayData(seu, assay = 'RNA', slot = 'counts')
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,cell_metadata = seu@meta.data, gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seu, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds)
cds <- learn_graph(cds,close_loop = F,use_partition = F) 

mytheme <-  theme(axis.title = element_text(size = 15, colour = "black"),
                  axis.text = element_text(size = 13, colour = "black"),
                  legend.text = element_text(size = 13, colour = "black"),
                  legend.title = element_text(size = 13, colour = "black"),
                  title = element_text(size = 15, colour = "black"))

p1 = plot_cells(cds, color_cells_by="sub_cell",label_groups_by_cluster = FALSE,
                group_label_size = 5, label_leaves = FALSE,label_branch_points = FALSE) + mytheme

p1
cds <- order_cells(cds)


p2 = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
                label_branch_points = FALSE) + mytheme
p2
patchwork::wrap_plots(p1, p2)
ggsave("B_sub_monocle3.pdf",width = 12,height = 5.5)

saveRDS(cds, "B_sub_monocle3.rds")
