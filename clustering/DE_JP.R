setwd("Master/bin/Merged_MAD/")
library(DESeq2)

load("model_JP.Rdata")
load("../../Data/Japanese data/exp_counts_gene_level_round.Rdata")

exp <- merge(exp, map_gene, by= "row.names", all.x=TRUE)
exp <- exp[!is.na(exp$Gene),]
exp[,c("Row.names")] <- NULL
### sum transcripts expression for each gene
exp <- aggregate(.~Gene, exp, sum)
row.names(exp) <- exp$Gene
exp[,1] <- NULL

groups <- as.data.frame(predict(model, what = "samples"))
colnames(groups) <- "cluster"
groups <- groups[order(match(row.names(groups), colnames(exp))),1, drop=FALSE]

exp <- as.data.frame(exp)

deseq <- DESeqDataSetFromMatrix(countData=exp, colData=groups, design = ~cluster)
deseq <- DESeq(deseq)
write.table(results(deseq, contrast=c("cluster","1","2")), file = "../../Results/Merged_MAD/Deseq/Cluster_1_2_DE_JP.txt", sep ="\t")
write.table(results(deseq, contrast=c("cluster","1","3")), file = "../../Results/Merged_MAD/Deseq/Cluster_1_3_DE_JP.txt", sep ="\t")
write.table(results(deseq, contrast=c("cluster","2","3")), file = "../../Results/Merged_MAD/Deseq/Cluster_2_3_DE_JP.txt", sep ="\t")

save(deseq, file = "deseq_TCGA.Rdata")
