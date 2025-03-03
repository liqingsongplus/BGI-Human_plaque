#!/usr/bin/Rscript
## Screening for effective SNP mutations

args<-commandArgs(T)
raw_file=as.character(args[1])  # processed.MAE_mito.rds 
anntation = args[2] # Two columns, separated by commas, with the first column beginning with: Cell classification _ Cell name.bam, eg:"../cellname.txt"

library(MultiAssayExperiment)
library(Matrix)
library("reshape2")
library("gplots")

"%ni%" <- Negate("%in%")
raw <- readRDS(raw_file)

covSE <- raw[["coverage"]]
allSE <- raw[["alleles"]]
baq <- assays(allSE)[["BAQ"]]
rm(raw)
covCell <-  Matrix::colMeans(assays(covSE)[["coverage"]])

# Allele Frequency
af <- assays(allSE)[["counts"]]/(assays(covSE)[["coverage"]][start(rowRanges(allSE)),] + 0.001)
rownames(af) <- paste0(data.frame(rowRanges(allSE))[,c(2)], "_", data.frame(rowRanges(allSE))[,c(7)])
cov <- assays(covSE)[["coverage"]][start(rowRanges(allSE)),]

master_df <- data.frame(cbind(data.matrix(af), data.matrix(cov), data.matrix(baq)))
master_df$position <- data.frame(rowRanges(allSE))[,c(2)]
master_df$altAllele <- data.frame(rowRanges(allSE))[,c(7)]
master_df$refAllele <- data.frame(rowRanges(allSE))[,c(6)]
pos=unique(master_df$position)


## quanlity control
num=dim(master_df)[2]/3-1
master_df2=master_df[,1:num]


myfun = function(x){
  a = as.numeric( gsub(" ", "", x[(num*2+1):(num*3)] ))
  b = as.numeric( gsub(" ", "", x[(num+1):(num*2)] ))
  if(median(melt(a)$value)<30 & median(melt(b)$value)<20 ){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

out = apply(master_df, 1, myfun)
master_df2 = master_df2[out,]

saveRDS(master_df2,file="master.rds")
