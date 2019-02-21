#setwd("C:/Users/kajetan.juszczak/Documents/Master/Results/Merged_MAD/Deseq/TCGA/")
setwd("/home/xnestea/Master/Results/Merged_MAD/Deseq/TCGA/")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("piano", version = "3.8")

library(piano)

l <- list.files(path = "./")
gset <- loadGSC("../../../../Data/h.all.v6.2.symbols.gmt")
for(i in l){
  setwd("/home/xnestea/Master/Results/Merged_MAD/Deseq/TCGA/")
  DE <- read.table(file = i)
  DE <- DE[complete.cases(DE),]
  name <- gsub("_DE.txt", "", i)
  DE=DE[ ,c('log2FoldChange','padj')]
  pval= as.matrix(DE[, 2]) #extract P as a matrix
  fc= as.matrix(DE[, 1])  #extract fold changes as a matrix
  row.names(pval)=row.names(DE)
  row.names(fc)=row.names(DE)
  gsaRes <- runGSA(pval,fc,gsc=gset, nPerm = 10000, adjMethod = "fdr")
  setwd("/home/xnestea/Master/Results/Merged_MAD/piano/TCGA/")
  GSAsummaryTable(gsaRes, save = TRUE, paste(name,"_piano_HM_TCGA.txt",sep=""))
}