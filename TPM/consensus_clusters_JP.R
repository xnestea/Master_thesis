### read in files related to one cluster
for(cluster in list_of_clusters){
  for (one_DE in all_possible){
    df <- read.table(paste(), sep = "\t")
    df <- df[c("padj" , "log2FoldChange")]
    df <- df[which(df$padj < 0.05),]
    temp <- rownames(df)
    savelist <- intersect(temp, savelist)
  write.table(savelist, file = paste(savelist))
  
}
}