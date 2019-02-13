setwd("Master/bin/")
load(file = "TPM/prep_exp_TCGA.Rdata")
load(file = "../Data/Japanese data/exp_tpm_gene_level.Rdata")
#symbol_exp<- as.data.frame(symbol_exp)
exp <- symbol_exp[which(rownames(symbol_exp) %in% rownames(exp)) ,]

save(exp, file = "TPM/prep_exp_JP.Rdata")
