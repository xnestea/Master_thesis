setwd("Master/bin/")

load("../Data/Japanese data/exp_tpm_gene_level.Rdata")
load("../Data/TCGA data/uni_transcript_tpm_exp_clinical_add1.Rdata")
load("TPM/map_gene.Rdata")

##########################################################################
#################################TCGA#####################################
##########################################################################
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
##########################################################################
####################################JP####################################
##########################################################################
exp_JP <- symbol_exp
MAD_JP <- apply(exp_JP, 1, mad)
exp_JP <- cbind(exp_JP, MAD_JP)
### order with respect to MAD
exp_JP <- exp_JP[order(-MAD_JP),]
exp_JP <- as.data.frame(exp_JP[1:1500, -ncol(exp_JP)])
##########################################################################
###################################MERGE##################################
##########################################################################
common <- rownames(exp)[rownames(exp) %in% rownames(exp_JP)]
##########################################################################
###############################DATA_PREP##################################
##########################################################################
exp <- exp[which(rownames(exp) %in% common),]
exp_JP <- exp_JP[which(rownames(exp_JP) %in% common),]
save(exp, file = "./Merged_MAD/Exp_TCGA.Rdata")
save(exp_JP, file = "./Merged_MAD/Exp_JP.Rdata")