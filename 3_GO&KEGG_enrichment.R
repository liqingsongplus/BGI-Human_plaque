## Enrichment analysis of a specified gene set
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)

outname = "Geneset"
geneset = readRDS("deg_marker.csv") # Read the differential gene set

## GO enrichment
go_Up <- enrichGO(geneset, OrgDb = "org.Hs.eg.db", ont="BP",keyType = 'SYMBOL',readable = T)

if(is.null(go_Up)){next}

go_Up = simplify(
  go_Up,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)

write.csv(go_Up, file = paste0(outname,'_GO.txt'))
dotplot(go_Up,showCategory=20,title = 'BP',label_format = 55) +  # orderBy= "",
  theme(axis.text.x = element_text(angle = 0))
ggsave( paste0(outname,"_GO.pdf"), width = 8, height = 8)

## KEGG enrichment
A = bitr(geneset, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL",'SYMBOL'), OrgDb="org.Hs.eg.db")[,2]
kk = enrichKEGG(A,keyType = "kegg",organism= "Hsa", qvalueCutoff = 0.05, pvalueCutoff=0.05)
kk <- setReadable(kk, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
hh <- as.data.frame(kk)
write.table(hh, paste0(outname,"_KEGG.txt"),sep = "\t",quote = F,row.names = F)

kk@result$Description = gsub(" -.*","", kk@result$Description)
dotplot(kk,showCategory=20,title = 'KEGG',label_format = 40, font.size = 10)
  theme(axis.text.x = element_text(angle = 0))
ggsave( paste0(outname,"_Kegg.pdf"), width = 6, height = 6)

