## The Stereo-seq data could be obtained from the CNGB Spatial Transcript Omics 
## DataBase (https://db.cngb.org/stomics/) with accession number: STT0000052.

# Set the spatial data to bin50

library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(ggplot2)
library(ggsci)


setwd("D:/work/3_4_human_plaque_2024/1_raw_sp")


A = dir("./1_raw_cellbin_gem/")  # gem file of stereo-seq

for(i in A){
  for(j in c(50)){
    
    infile <- i
    bs <- j  # 50
    
    pro = gsub("_.*","", infile)
    
    dat <- fread(file = paste0("./1_raw_cellbin_gem//", infile))
    dat$label = NULL
    dat$tag = NULL
    
    if(length(grep("MIDCounts|MIDCount",colnames(dat))>0)){
      colnames(dat) <- gsub("MIDCounts|MIDCount","UMICount",colnames(dat))}
    out <- as.data.frame(dat)
    
    # 标准化芯片距离，获得bin
    dat$x <- trunc((dat$x - min(dat$x)) / bs + 1)
    dat$y <- trunc((dat$y - min(dat$y)) / bs + 1)
    
    out <- cbind(dat$y,dat$x,out)
    colnames(out)[1:2] <- c(paste0("bin",bs,".y"),paste0("bin",bs,".x"))
    # fwrite(out,paste0(pro,"_bin",bs,"_information.txt"),col.names=T,row.names=F,sep="\t",quote=F)
    
    dat <- dat[, sum(UMICount), by = .(geneID, x, y)] 
    dat$bin_ID <- max(dat$x) * (dat$y - 1) + dat$x 
    bin.coor <- dat[, sum(V1), by = .(x, y)]   
    

    geneID <- seq(length(unique(dat$geneID)))  
    hash.G <- data.frame(row.names = unique(dat$geneID), values = geneID) 
    gen <- hash.G[dat$geneID, 'values'] 
    

    bin_ID <- unique(dat$bin_ID) 
    hash.B <- data.frame(row.names = sprintf('%d', bin_ID), values = bin_ID)
    bin <- hash.B[sprintf('%d', dat$bin_ID), 'values']
    
    cnt <- dat$V1
    

    ##
    tissue_lowres_image <- matrix(1, max(bin.coor$y), max(bin.coor$x))  
    
    tissue_positions_list <- data.frame(row.names = paste('BIN', rownames(hash.B), sep = '.'),
                                        tissue = 1, 
                                        row = bin.coor$y, col = bin.coor$x,
                                        imagerow = bin.coor$y, imagecol = bin.coor$x)
    
    scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1,
                                     tissue_hires_scalef = 1,
                                     tissue_lowres_scalef = 1))  
    
 
    mat <- sparseMatrix(i = gen, j = bin, x = cnt)  # i为基因，j为细胞，x为（i,j）处表达量
    
    rownames(mat) <- rownames(hash.G)
    colnames(mat) <- paste('BIN', sprintf('%d', seq(max(hash.B[, 'values']))), sep = '.')
    
    ############################## 02. creat Spatial Object  ##############################
    seu <- CreateSeuratObject(mat, project = 'Spatial', assay = 'Spatial',min.cells=3, min.features=10)
    seu = subset(seu, nFeature_Spatial > 20)
    
    generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) {
      if (filter.matrix) {
        tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
      }
      
      unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
      
      spot.radius <- unnormalized.radius / max(dim(image))
      
      return(new(Class = 'VisiumV1', 
                 image = image, 
                 scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
                                              fiducial = scale.factors$fiducial_diameter_fullres, 
                                              hires = scale.factors$tissue_hires_scalef, 
                                              lowres = scale.factors$tissue_lowres_scalef), 
                 coordinates = tissue.positions, 
                 spot.radius = spot.radius))
    }
    
    spatialObj <- generate_spatialObj(image = tissue_lowres_image, 
                                      scale.factors = fromJSON(scalefactors_json), 
                                      tissue.positions = tissue_positions_list)
    
    spatialObj <- spatialObj[Cells(seu)]
    DefaultAssay(spatialObj) <- 'Spatial'
    
    seu[['slice1']] <- spatialObj
    
    coor = GetTissueCoordinates(object = seu[["slice1"]])
    seu = AddMetaData(seu, metadata = coor)
    ggplot(seu@meta.data, aes(x = imagerow, y = imagecol, color = nFeature_Spatial)) + geom_point(size = 0.01)
    

    ##############################  03. Spatial Analyse  ##############################
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")  # 每个细胞检测到的线粒体基因比例
    
    Q1 <- quantile(seu$nFeature_Spatial)[2]
    Q3 <- quantile(seu$nFeature_Spatial)[4]
    upper <- as.numeric(Q3+1.5*(Q3-Q1))
    lower <- as.numeric(Q1-1.5*(Q3-Q1))
    
    p1 <- VlnPlot(seu, features=c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol=3, pt.size=0)+theme(axis.text.x=element_text(angle=0,size=0),axis.title.x=element_text(angle=20,size=8))
    print(p1)
    p2 <- ggplot(seu@meta.data,aes(x=nFeature_Spatial)) +geom_density(colour="black") + theme_classic() + 
      theme(plot.title=element_text(hjust=0.5,size=18, face="bold.italic"), legend.position="none",axis.title=element_text(size=15, face="bold.italic"),axis.text.x=element_text(size=12),axis.ticks.x=element_blank()) + 
      geom_vline(aes(xintercept=100,colour="#999999",linetype="twodash")) + geom_vline(aes(xintercept=200,colour="#999999",linetype="twodash"))+
      geom_vline(aes(xintercept=300,colour="#999999",linetype="twodash"))+geom_vline(aes(xintercept=500,colour="#999999",linetype="twodash")) + 
      xlim(min(seu@meta.data$nFeature_Spatial),max(seu@meta.data$nFeature_Spatial))
    print(p2)
    plot1 <- FeatureScatter(seu, feature1 = "nCount_Spatial", feature2 = "percent.mt") + NoLegend()
    plot2 <- FeatureScatter(seu, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") + NoLegend()
    
    pdf(file=paste0(pro, "_QC.pdf"), width=20,height=20)
    cowplot::plot_grid(p1, p2,plot1,plot2, rel_widths=c(1,1,1,1), ncol=2)
    dev.off()
    
    seu$orig.ident = pro
    
    saveRDS(seu,file=paste0(pro,"_bin",bs,"_preQC.rds"))
    
    rm(list = ls()[! ls() %in% c("A", "i")])
    gc()
  }
  
}
