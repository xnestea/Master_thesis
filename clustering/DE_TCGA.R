setwd("/home/xnestea/Master/bin/Merged_MAD/")

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

install.packages("NMF")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")

library(NMF)
library(DESeq2)

load("../TPM/map_gene.Rdata")
load("model_TCGA.Rdata")
load("../../Data/TCGA data/uni_transcript_counts_exp_clinical_add1.Rdata")

exp <- merge(exp, map_gene, by= "row.names", all.x=TRUE)
exp <- exp[!is.na(exp$Gene),]
exp[,c("Row.names")] <- NULL
### sum transcripts expression for each gene
exp <- aggregate(.~Gene, exp, sum)
row.names(exp) <- exp$Gene
exp[,1] <- NULL
#save_exp <- exp
#exp <- save_exp
### retrieve clusters
groups <- as.data.frame(predict(model, what = "samples"))
colnames(groups) <- "cluster"
groups <- groups[order(match(row.names(groups), colnames(exp))),1, drop=FALSE]

exp[] <- sapply(exp, as.integer)

deseq <- DESeqDataSetFromMatrix(countData=exp, colData=groups, design = ~cluster)
deseq <- DESeq(deseq)
write.table(results(deseq, contrast=c("cluster","1","2")), file = "../../Results/Merged_MAD/Deseq/TCGA/Cluster_1_2_DE.txt", sep ="\t")
write.table(results(deseq, contrast=c("cluster","1","3")), file = "../../Results/Merged_MAD/Deseq/TCGA/Cluster_1_3_DE.txt", sep ="\t")
write.table(results(deseq, contrast=c("cluster","2","3")), file = "../../Results/Merged_MAD/Deseq/TCGA/Cluster_2_3_DE.txt", sep ="\t")

save(deseq, file = "deseq_TCGA.Rdata")