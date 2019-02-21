setwd("C:/Users/kajetan.juszczak/Documents/Master/bin/")
load("./Merged_MAD/model_JP.Rdata")

JP_clus <- read.csv("../Data/ng.2699-S2.csv")
rownames(JP_clus) <- JP_clus$Supplementary.table.1..Characteristics.of.the.patinets
JP_clus <- JP_clus[which(JP_clus$X.11 == "ccA" | JP_clus$X.11 == "ccB" ),]

groups <- as.data.frame(predict(model, what = "samples"))
colnames(groups) <- "cluster"
r <- str_locate(rownames(groups), pattern = "cc.*-")
r <- as.data.frame(r)
rownames(groups) <- substr(x = rownames(groups),r$start,r$end - 1)

df <- merge(JP_clus, groups, by = "row.names")
rownames(df) <- df$Row.names
df <- df[,c("cluster", 'X.11')]
df$count <- 1
agg <- aggregate(count ~ ., df, FUN = sum)
agg <- as.data.frame(agg[order(agg$cluster),])

write.table(agg, file = "../Results/Merged_MAD/DB_sub_comp_JP.txt", sep = "\t")
