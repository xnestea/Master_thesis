library(stringr)
library(NMF)
library(plyr)


setwd("c:/Users/kajetan.juszczak/Documents/Master/")
### load NMF model
load("./bin/Merged_MAD/model_JP.Rdata")
#format ensebl - symbol dicitonary

groups <- as.data.frame(predict(model, what = "samples"))
colnames(groups) <- "cluster"
groupsJ <- groups[order(match(row.names(groups), colnames(exp))),1, drop=FALSE]

load("./bin/Merged_MAD/model_TCGA.Rdata")
#format ensebl - symbol dicitonary
groups <- as.data.frame(predict(model, what = "samples"))
colnames(groups) <- "cluster"
groupsT <- groups[order(match(row.names(groups), colnames(exp))),1, drop=FALSE]

count(groupsJ, "cluster")
count(groupsT, "cluster")