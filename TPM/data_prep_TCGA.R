setwd("Master/bin/")

library(SummarizedExperiment)
library(NMF)
library(biomaRt)

load("../../Data/TCGA data/uni_transcript_tpm_exp_clinical_add1.Rdata")

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
### for all ensemble genes in dataframe add information about transcript biotype and gene symbol
map_gene = getBM(attributes=c('ensembl_transcript_id','ensembl_gene_id',"transcript_biotype",'external_gene_name'), 
                 filters = 'ensembl_transcript_id', values=rownames(exp), mart=ensembl)
### select only protein coding.
#save_mp <- map_gene
#map_gene <- save_mp
map_gene = map_gene[which(map_gene$transcript_biotype=="protein_coding"),]
### change name of column
colnames(map_gene)[colnames(map_gene)=="external_gene_name"] <- "Gene"
### set correct rownames
rownames(map_gene) <- map_gene[,1]
### remove redundand columns
map_gene[,c(1:3)] <- NULL
### not used
#test = save_exp[which(save_exp$Gene == "A1BG"),]
#exp2<- exp[which(rownames(exp) %in% test),]
#save_exp <- exp
#exp <- save_exp
### change ensembl to gene symbol
exp <- merge(exp, map_gene, by= "row.names", all.x=TRUE)
exp <- exp[!is.na(exp$Gene),]
exp[,c("Row.names")] <- NULL
exp <- aggregate(.~Gene, exp, sum)
#save_agg_exp <- exp
#exp <- save_agg_exp
row.names(exp) <- exp$Gene
exp[,1] <- NULL
### calculate MAD
MAD <- apply(exp, 1, mad)
exp <- cbind(exp, MAD)
### order with respect to MAD
exp <- exp[order(-MAD),]
### take top 1500
exp <- exp[1:1500, -ncol(exp)]

save(exp, file = "TPM/prep_exp_TCGA.Rdata")
