setwd("C:/Users/kajetan.juszczak/Documents/Master/Results/Merged_MAD/Deseq/")

t <- list.files(path = "TCGA/")
j <- list.files(path = "JP/")

for(i in t){
  setwd("C:/Users/kajetan.juszczak/Documents/Master/Results/Merged_MAD/Deseq/TCGA/")
  temp_t <- read.table(file = i)
  temp_t <- temp_t[order(temp_t$padj),]
  sum_gene <- length(rownames(temp_t))
  for(x in j){
    setwd("C:/Users/kajetan.juszczak/Documents/Master/Results/Merged_MAD/Deseq/JP/")
    temp_j <- read.table(file = x)
    temp_j <- temp_j[order(temp_j$padj),]
    T1 <- temp_t[1:2000,]
    J1 <- temp_j[1:2000,]
    a <- rownames(T1)
    b <- rownames(J1)
    intersect <- intersect(a, b)
    ratio <- length(intersect) / 2000
    T1 <- T1[intersect,]
    J1 <- J1[intersect,]
    s <- sign(T1$log2FoldChange) == sign(J1$log2FoldChange)
    same <- length(s[s == T])
    v <- same / length(intersect)
    #print(s)
    

    print(paste(i,x,ratio, v))
    T1 <- T1[intersect,]
    J1 <- J1[intersect,]
    s <- sign(T1$log2FoldChange) == sign(J1$log2FoldChange)
    #print(s)
    same <- length(s[s == T])
    p <- same / length(intersect)
    if(p < 0.05 | p > 0.95){
      res <- 1 - phyper(length(intersect) - 1, 2000, sum_gene - 2000, 2000)
      #print(res)
      #p_hyper=1-phyper(dim(over_up)[1]+dim(over_down)[1]-1,num_deg_1,num_deg_2)
    }
  }
}
