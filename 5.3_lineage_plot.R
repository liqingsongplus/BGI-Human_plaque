## Mapping based on mitochondrial snp mutations

library(Seurat)
library(gplots)
library("ggplot2")
library("Rtsne")
library(ape)
library(tidytree)
library(ggtree)
library(tidyr)

setwd("/B_puxi/")  # master.rds file path
ll = dir("./")

for(rds in ll){
  df = readRDS(rds)
  i = gsub("_master.rds","", rds)
  colnames(df) = gsub("cell_","",colnames(df))
  sds=apply(df,1,sd)
  df=df[sds>0.3,]
  
  ## TSNE
  t=Rtsne(t(df),perplexity=5,check_duplicates =F, normalize = F)
  ty=as.data.frame(t$Y)
  rownames(ty) = colnames(df)
  names(ty)=c("TSNE_1","TSNE_2")
  write.csv(ty, paste0(i,"_tsne.csv"))
  # ty = read.csv(paste0(i,"_tsne.csv"),row.names = 1)
  ty$group = gsub("_.*","",rownames(ty))
  ggplot(data=ty,aes(x=`TSNE_1`,y=`TSNE_2`,color=group))+
    geom_point()+theme_bw()+
    theme(panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line=element_line(color="black"))
  ggsave(paste0(i,"_tsne_sd_0.3.pdf"),width = 7,height = 6)

  ## Evolutionary tree
  hc = hclust(dist(t(df)))
  groupinfo = split(rownames(ty), ty$group) 
  tree = tidytree::groupOTU(as.phylo(hc), groupinfo)
  
  ggtree(tree,layout = "circular", aes(color = group), ladderize = F, branch.length = "none",size = 0.05) +
    labs(color = NULL) 
  ggsave(paste0(i,"_tree.pdf"), width = 7,height = 7)
  saveRDS(tree, paste0(i,"_tree.rds"))
  
  rm(tree, df, group, groupinfo, hc, cn, sds)
}
