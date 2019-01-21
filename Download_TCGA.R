if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

ccRCC <- GDCquery(project = "TCGA-KIRC",
                  data.category = "Transcriptome Profiling",
                  legacy = FALSE,
                  workflow.type = "HTSeq - Counts",
                  sample.type = c("Solid Tissue Normal"))

setwd("Master/Scripts/")
GDCdownload(ccRCC)