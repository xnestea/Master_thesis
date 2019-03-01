setwd("C:/Users/kajetan.juszczak/Documents/Master/bin/")
load("./Merged_MAD/model_TCGA.Rdata")

TCGA_clus <- read.csv("../Data/Data_file_S9_mRNA_miRNA_cluster_assignments.xlsx")
rownames(TCGA_clus) <- TCGA_clus$Patient

groups <- as.data.frame(predict(model, what = "samples"))
colnames(groups) <- "cluster"

df <- merge(TCGA_clus, groups, by = "row.names")
rownames(df) <- df$Patient
df <- df[,c("cluster", "mRNA_cluster")]
df$count <- 1
agg <- aggregate(count ~ ., df, FUN = sum)
agg <- as.data.frame(agg[order(agg$cluster),])

write.table(agg, file = "../Results/Merged_MAD/DB_sub_comp_TCGA.txt", sep = "\t")
