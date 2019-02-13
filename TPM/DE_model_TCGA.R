setwd("Master/bin/")

library(DESeq2)

load("../../bin/TPM/map_gene.Rdata")
load("../../bin/TPM/model_TCGA.Rdata")
load("../../Data/TCGA data/uni_transcript_counts_exp_clinical_add1.Rdata")
### set gene names as rownames
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
res <- results(deseq)
write.table(res, file = "../Results/TPM/DE_TCGA.txt", sep = "\t")