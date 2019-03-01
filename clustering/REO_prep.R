library(stringr)
library(NMF)

setwd("c:/Users/kajetan.juszczak/Documents/Master/")
### load exp matrix
load("./Data/Japanese data/exp_tpm_gene_level.Rdata")
### load NMF model
load("./bin/Merged_MAD/model_JP.Rdata")
#format ensebl - symbol dicitonary
rownames(ensem_of_symbol_exp) <- ensem_of_symbol_exp[,2] 
ensem_of_symbol_exp<- ensem_of_symbol_exp[, 1]
### merge to convert expression from gene_symbols to ensembl id
exp <- merge(symbol_exp, as.data.frame(ensem_of_symbol_exp), by=0, all=TRUE)
### set ensembl IDs as rownames
rownames(exp) <- exp$ensem_of_symbol_exp
### remove redundant rows
exp[, c("Row.names", "ensem_of_symbol_exp")] <- NULL
### cut begining of ensembl IDs - necessary for REO
rownames(exp)<- substr(rownames(exp),5, 25)

### retrive clusters based on model and order them to be the same as exp matrix
groups <- as.data.frame(predict(model, what = "samples"))
### correct name of the column
colnames(groups) <- "cluster"
### make sure exp and group samples are in the same order
groups <- groups[order(match(row.names(groups), colnames(exp))),1, drop=FALSE]
### for REO in matlab rownames must be numeric so I cut non numeric part
#rownames(groups)<- substr(rownames(groups),5, 15)
#groups = as.matrix(groups)
#write.table(groups, file = "Data/cluster_table.txt", sep = "\t", col.names = F, row.names = F, quote = F)

### select expression data for each cluster and save it to separate files
for(i in list("1", "2", "3")){
  ### find names of samples in current cluster
  cluster_samples <- rownames(groups)[which(groups$cluster == i)]
  ### filter expresion table to samples from current cluster
  cexp = exp[, which(colnames(exp) %in% cluster_samples)]
  ### write it to file
  write.table(cexp, file = paste0("./Data/TPM_JP_", i,".txt"), sep = "\t")
}


library(biomaRt)
load("./Data/TCGA data/uni_transcript_tpm_exp_clinical_add1.Rdata")
load("./bin/Merged_MAD/model_TCGA.Rdata")

### download ensembl query
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
### for all ensemble genes in dataframe add information about transcript biotype
map_gene = getBM(attributes=c('ensembl_transcript_id','ensembl_gene_id',"transcript_biotype"), 
                 filters = 'ensembl_transcript_id', values=rownames(exp), mart=ensembl)
### backup
#save_mp <- map_gene
#map_gene <- save_mp
### select only protein coding.
map_gene = map_gene[which(map_gene$transcript_biotype=="protein_coding"),]
### change name of column
colnames(map_gene)[colnames(map_gene)=="ensembl_gene_id"] <- "Gene"
### set correct rownames
rownames(map_gene) <- map_gene[,1]
### remove redundand columns
map_gene[,c(1,3)] <- NULL

### add gene column to exp dataframe
exp <- merge(exp, map_gene, by= "row.names", all.x=T)
#exp2 <- exp
#exp <- exp2
### drop rows without gene symbol
exp <- exp[!is.na(exp$Gene),]
### Remove redundant column
exp[,c("Row.names")] <- NULL
### sum transcripts expression for each gene
exp <- aggregate(.~Gene, exp, sum)
### Set Gene ID as rownames
row.names(exp) <- exp$Gene
### Remove ths column after
exp[,1] <- NULL
### cut begining of ensembl IDs - necessary for REO
rownames(exp)<- substr(rownames(exp),5, 25)

### retrive clusters based on model and order them to be the same as exp matrix as above
groups <- as.data.frame(predict(model, what = "samples"))
colnames(groups) <- "cluster"
groups <- groups[order(match(row.names(groups), colnames(exp))),1, drop=FALSE]
#groups = as.matrix(groups)
#write.table(groups, file = "Data/cluster_table_TCGA.txt", sep = "\t", col.names = F, row.names = F, quote = F)

### select expression data for each cluster and save it to separate files as above
for(i in list("1", "2", "3")){
  cluster_samples <- rownames(groups)[which(groups$cluster == i)]
  cexp = exp[, which(colnames(exp) %in% cluster_samples)]
  write.table(cexp, file = paste0("./Data/TPM_TCGA_", i,".txt"), sep = "\t")
}

