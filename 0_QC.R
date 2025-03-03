## Quality control of single cell data in 10× format.
## The integrated single-cell data can be downloaded from CNSA using accession number: CNP0004894. (https://db.cngb.org/cnsa)


library(ggplot2)
library(Seurat)
library(DoubletFinder)

setwd("D:/work/3_4_human_plaque_2024/11_norm_PLA/")  # Set the location of single cell data for different samples.
dirlist = dir()
for(i in dirlist){
  setwd(paste0("D:/work/3_4_human_plaque_2024/11_norm_PLA/",i))
  seu = Read10X("SpliceMatrix/", gene.column = 1)  # 每个样本的10×格式单细胞文件储存位置
  seu <- CreateSeuratObject(counts = seu, project = "RNA", min.cells = 3, min.features = 10)
  seu$orig.ident = i
  
  seu = RenameCells(seu,new.names=paste(i,colnames(seu),sep="_"))
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  dir.create("figures_QC")
  pdf(paste0("figures_QC/preQC.pdf"),width=10,height=8)
  print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
  print(ggplot(seu@meta.data,aes(x=nFeature_RNA)) +geom_density(colour="black") + theme_classic() + 
          theme(plot.title=element_text(hjust=0.5,size=18, face="bold.italic"), legend.position="none",axis.title=element_text(size=15, face="bold.italic"),axis.text.x=element_text(size=12),axis.ticks.x=element_blank()) + 
          geom_vline(xintercept=c(100,200,300,500,mean(seu$nFeature_RNA)),color="red",linetype="twodash") + 
          xlim(min(seu@meta.data$nFeature_RNA),max(seu@meta.data$nFeature_RNA)))
  
  plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  dev.off()

  ###### Define Find_doublet function ########
  Find_doublet <- function(data){
    set.seed(123)
    sweep.res.list <- paramSweep_v3(data, PCs = 1:20, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    doublets.percentage = 0.05
    nExp_poi <- round(as.numeric(doublets.percentage)*ncol(data))
    p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
    data <- doubletFinder_v3(data, PCs = 1:20, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
    # data<-subset(data,subset=doublet_info=="Singlet")
    return(data)
  }
  
  seu = seu |> NormalizeData() |>
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) |>
    ScaleData() |>
    RunPCA() |>
    RunUMAP(dims = 1:30) |>
    Find_doublet()
  
  write.table(seu@meta.data,paste0("doublets_info.txt"),sep="\t",quote=FALSE)
  seu <- subset(seu,subset=doublet_info=="Singlet" & percent.mt<20 & nFeature_RNA>200)
  
  pdf(paste0("figures_QC/QC.pdf"),width=10,height=8)
  print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
  print(ggplot(seu@meta.data,aes(x=nFeature_RNA)) +geom_density(colour="black") + theme_classic() + 
          theme(plot.title=element_text(hjust=0.5,size=18, face="bold.italic"), legend.position="none",axis.title=element_text(size=15, face="bold.italic"),axis.text.x=element_text(size=12),axis.ticks.x=element_blank()) + 
          geom_vline(xintercept=c(100,200,300,500,700,1000),color="red",linetype="twodash") + 
          xlim(min(seu@meta.data$nFeature_RNA),max(seu@meta.data$nFeature_RNA)))
  
  plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  dev.off()
  saveRDS(seu, paste0(i,"_QC.rds"))   # The quality control results for each single cell sample are stored in sample name +_QC.rds format
  rm(seu, plot1, plot2)
  gc()
}