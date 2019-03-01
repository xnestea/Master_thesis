setwd("C:/Users/kajetan.juszczak/Documents/Master/Results/Merged_MAD/Deseq/")



##### AIM OF THE SCRIPT IS TO FIND RESPECTIVE CLUSTERS BETWEEN 2 DATABESES
##### AND TO SEE IF SUCH PAIRS ARE RELIABLE 



### list all deseq files from folders
### for JP and TCGA data respectively
t <- list.files(path = "TCGA/")
j <- list.files(path = "JP/")

### for each deseq file in TCGA
for(i in t){
  ### set corect wd
  setwd("C:/Users/kajetan.juszczak/Documents/Master/Results/Merged_MAD/Deseq/TCGA/")
  ### read in file and order its entries with repsect to p.adj - increasing order
  temp_t <- read.table(file = i)
  temp_t <- temp_t[order(temp_t$padj),]
  ### for later tests save number of genes in total
  sum_gene <- length(rownames(temp_t))
  ### for every deseq in JP data - ensure comparison each vs each
  for(x in j){
    ### set corect wd
    setwd("C:/Users/kajetan.juszczak/Documents/Master/Results/Merged_MAD/Deseq/JP/")
    ### read in file and order its entries with repsect to p.adj - increasing order
    temp_j <- read.table(file = x)
    temp_j <- temp_j[order(temp_j$padj),]
    ### select top 2000 genes with respect to p.adj value
    T1 <- temp_t[1:2000,]
    J1 <- temp_j[1:2000,]
    ### save names of genes as a and b variable
    a <- rownames(T1)
    b <- rownames(J1)
    ### list of comon gene names
    intersect <- intersect(a, b)
    ### how many genes are overlaping with respect to the total number of selected genes = 2000 
    ratio <- length(intersect) / 2000
    ### Table for each 2000 top gene tables for respective database
    ### select only genes that are overlaping
    T1 <- T1[intersect,]
    J1 <- J1[intersect,]
    ### calculate the percentage of overlapping genes having the same sign
    s <- sign(T1$log2FoldChange) == sign(J1$log2FoldChange)
    same <- length(s[s == T])
    v <- same / length(intersect)
    #print(s)
    ### print to see which deseqs are reversed
    print(paste(i,x,ratio, v))
    T1 <- T1[intersect,]
    J1 <- J1[intersect,]
    s <- sign(T1$log2FoldChange) == sign(J1$log2FoldChange)
    #print(s)
    same <- length(s[s == T])
    p <- same / length(intersect)
    ### for comparisons with sufficient concordance in directions of expression
    if(p < 0.05 | p > 0.95){
      ### calculate hypergeometric test result
      ### p-value = chance of result being random
      res <- 1 - phyper(length(intersect) - 1, 2000, sum_gene - 2000, 2000)
      #print(res)
      #p_hyper=1-phyper(dim(over_up)[1]+dim(over_down)[1]-1,num_deg_1,num_deg_2)
    }
  }
}
