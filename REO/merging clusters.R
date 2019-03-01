library(stringr)
setwd("C:/Users/kajetan.juszczak/Documents/Master/")

### read in expression tables for both databases
TCGA <- read.table("Data/TPM_TCGA.txt", sep = "\t")
JP <- read.table("Data/TPM_JP.txt", sep = "\t")
### merge into one table 
mrg <- merge(TCGA, JP, by = 0)
## 75% of the sample size
smp_size <- floor(0.75 * ncol(mrg))

## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(ncol(mrg)), size = smp_size)

train <- mrg[, train_ind]
test <- mrg[, -train_ind]
### correct rownames
rownames(mrg) = mrg[,1]
mrg[,1] <- NULL

write.table(mrg, file = "Data/merged_exp_train.txt", sep = "\t")
### load 2 models and get clusters data
load("./bin/Merged_MAD/model_JP.Rdata")
model_JP <- model
groups <- as.data.frame(predict(model_JP, what = "samples"))
colnames(groups) <- "cluster"
groups_JP <- groups[order(match(row.names(groups), colnames(exp))),1, drop=FALSE]
load("./bin/Merged_MAD/model_TCGA.Rdata")

groups <- as.data.frame(predict(model, what = "samples"))
colnames(groups) <- "cluster"
groups <- groups[order(match(row.names(groups), colnames(exp))),1, drop=FALSE]

### translate JP clusters to TCGA clusters
groups_JP$new <- 0
groups_JP[which(groups_JP$cluster == 1),]$new <- 2
groups_JP[which(groups_JP$cluster == 2),]$new <- 3
groups_JP[which(groups_JP$cluster == 3),]$new <- 1
groups_JP$cluster <- groups_JP$new
groups_JP$new <- NULL

### megre cluster table
clusters <- rbind(groups, groups_JP)
rownames(clusters) <- gsub('-', '.', rownames(clusters))
### write as table
write.table(clusters, file = "Data/merged_clusters.txt",sep = "\t")

for(i in list("1", "2", "3")){
  ### find names of samples in current cluster
  cluster_samples <- rownames(clusters)[which(clusters$cluster == i)]
  ### filter expresion table to samples from current cluster
  cexp = mrg[, which(colnames(mrg) %in% cluster_samples)]
  ### write it to file
  write.table(cexp, file = paste0("./Data/TPM_Merged_train", i,".txt"), sep = "\t")
}
