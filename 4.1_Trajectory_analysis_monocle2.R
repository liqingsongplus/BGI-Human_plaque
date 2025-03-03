## A pseudo-time sequence analysis of a given cell type was performed using monocle2

library(monocle)
library(Seurat)
library(tidyverse)

seu = readRDS('B_sub.rds')  # Import the single-cell RNA sequencing data stored in Seurat object format to perform trajectory analysis.

data <- GetAssayData(seu, assay = "RNA",slot = 'counts')
pd <- new('AnnotatedDataFrame', data = seu@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)

cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds) # normalize for differences in mRNA recovered across cells
cds <- estimateDispersions(cds) # the mean-variance model learning by estimateDispersions()

diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~sub_cell",cores = 20) # Differential genes were obtained by grouping them according to the specified cell type
saveRDS(diff_test_res, "diff_gene.rds")
write.table(diff_test_res,"diff_gene.txt",sep = "\t",quote = F,row.names = F)

ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))

cds <- setOrderingFilter(cds, ordering_genes)

plot_ordering_genes(cds)
ggsave("plot_ordering_genes.png",width = 7,height = 6)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

cds <- orderCells(cds)
saveRDS(cds,paste0("orderCells.rds"))

mytheme <-  theme(axis.title = element_text(size = 15, colour = "black"),
                  axis.text = element_text(size = 13, colour = "black"),
                  legend.text = element_text(size = 13, colour = "black"),
                  legend.title = element_text(size = 13, colour = "black"),
                  strip.text = element_text(size=13),
                  title = element_text(size = 15, colour = "black")) 

mygroup = "sub_cell" # set group
plot_cell_trajectory(cds, color_by = mygroup) + mytheme + labs(colour = mygroup) 
ggsave('cell_traject.pdf',width = 6,height = 6)
plot_cell_trajectory(cds, color_by = mygroup,show_state_number = F,cell_size = 0.1 ) + facet_wrap(mygroup, ncol = 4) + mytheme  + labs(colour = "Cell type")
ggsave('cell_split.pdf',width = 12,height = 10)
plot_cell_trajectory(cds, color_by = "Pseudotime") + mytheme + scale_color_viridis_c()
ggsave('Pseudotime.pdf',width = 6,height = 6)
plot_cell_trajectory(cds, color_by = "State") + mytheme
ggsave('State.pdf',width = 6,height = 6)