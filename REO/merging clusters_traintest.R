library(stringr)
setwd("C:/Users/kajetan.juszczak/Documents/Master/")

### read in expression tables for both databases
TCGA <- read.table("Data/TPM_TCGA.txt", sep = "\t")
JP <- read.table("Data/TPM_JP.txt", sep = "\t")
### merge into one table 
mrg <- merge(TCGA, JP, by = 0)
### correct rownames
rownames(mrg) = mrg[,1]
mrg[,1] <- NULL
#################TRAINTEST_SPLIT############################
## 75% of the sample size
smp_size <- floor(0.70 * ncol(mrg))

## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(ncol(mrg)), size = smp_size)

train <- mrg[, train_ind]
test <- mrg[, -train_ind]

write.table(train, file = "Data/merged_exp_train.txt", sep = "\t")
write.table(test, file = "Data/merged_exp_test.txt", sep = "\t")

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
train_c <- as.data.frame(clusters[which(rownames(clusters) %in% colnames(train)), ,FALSE])
test_c <- as.data.frame(clusters[which(rownames(clusters) %in% colnames(test)), ,FALSE])
write.table(train_c, file = "Data/merged_clusters_train.txt",sep = "\t")
write.table(test_c, file = "Data/merged_clusters_test.txt",sep = "\t")
################TRAIN###################################################
for(i in list("1", "2", "3")){
  ### find names of samples in current cluster
  cluster_samples <- rownames(train_c)[which(train_c$cluster == i)]
  ### filter expresion table to samples from current cluster
  cexp = train[, which(colnames(train) %in% cluster_samples)]
  ### write it to file
  write.table(cexp, file = paste0("./Data/TPM_Merged_train_", i,".txt"), sep = "\t")
}
###############TEST####################################################
for(i in list("1", "2", "3")){
  ### find names of samples in current cluster
  cluster_samples <- rownames(test_c)[which(test_c$cluster == i)]
  ### filter expresion table to samples from current cluster
  cexp = test[, which(colnames(test) %in% cluster_samples)]
  ### write it to file
  write.table(cexp, file = paste0("./Data/TPM_Merged_test_", i,".txt"), sep = "\t")
}
