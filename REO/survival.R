library(NMF)
library(stringr)
library(survival)
setwd("C:/Users/kajetan.juszczak/Documents/Master/")

###################FOR MERGED DATA########################
### read cluster - sample info
###cluster_data <- read.table("./Data/merged_clusters.txt")

### format to match with survival info
### TIP cant use "." because it means any character in regexp I guess
###rownames(cluster_data) = gsub('[.]', '-', rownames(cluster_data))

#####################TCGA DATA######################
### load cluster anotation
load("bin/Merged_MAD/model_TCGA.Rdata")
TCGA_clusters <- as.data.frame(predict(model, what = "samples"))
colnames(TCGA_clusters) <- "cluster"
### MAIN AIM HERE:select only columns I need from clinical_exp
load("Data/TCGA data/uni_transcript_counts_exp_clinical_add1.Rdata")
### to index properly
TCGA_survival <- as.data.frame(clinical_exp)
### set barcodes as rownames
rownames(TCGA_survival) <- TCGA_survival$bcr_patient_barcode
### select interesting columns
TCGA_survival <- TCGA_survival[,c("days_to_death", "days_to_last_follow_up", "vital_status")]
### stupidest thing in the word. if converting factors have to use as numeric + as character otherwise it wont work
TCGA_survival$days <- ifelse(is.na(TCGA_survival$days_to_death), 
                            as.numeric(as.character(TCGA_survival$days_to_last_follow_up)), 
                            as.numeric(as.character(TCGA_survival$days_to_death)))
### need to be converted to 0 and 1 to work for survfit
TCGA_survival$vital_status <- ifelse(TCGA_survival$vital_status == "dead",1,0)
TCGA <- merge(TCGA_survival, TCGA_clusters, by = 0)
rownames(TCGA) <- TCGA$Row.names
TCGA[,1:3] <- NULL
####################JP DATA#########################
### load cluster anotation
load("bin/Merged_MAD/model_JP.Rdata")
JP_clusters <- as.data.frame(predict(model, what = "samples"))
colnames(JP_clusters) <- "cluster"
### cluster numbers to match TCGA
JP_clusters$new <- 0
JP_clusters[which(JP_clusters$cluster == 1),]$new <- 2
JP_clusters[which(JP_clusters$cluster == 2),]$new <- 3
JP_clusters[which(JP_clusters$cluster == 3),]$new <- 1
JP_clusters$cluster <- JP_clusters$new
JP_clusters$new <- NULL
### not sure if can be done in 1 line - format to match with sample names from survival sheet
rownames(JP_clusters) <- str_extract(rownames(JP_clusters), regex("(c.*-)"))
rownames(JP_clusters) <- substr(rownames(JP_clusters), 1, nchar(rownames(JP_clusters))-1)
### get survival data - sample relationship for JP
JP_survival <- read.csv(file = "Data/ng.2699-S2.csv", header = T)
### select only columns I need
rownames(JP_survival) = JP_survival$sample.ID
JP_survival <- JP_survival[,c("outcome","observation.period..month.")]
### merge both
JP <- merge(JP_survival, JP_clusters, all.y = T, by = 0)
rownames(JP) <- JP$Row.names
JP$Row.names <- NULL
colnames(JP) <- colnames(TCGA)
JP$days <- JP$days *30
### need to be converted to 0 and 1 to work for survfit
JP$vital_status <- ifelse(JP$vital_status == "dead",1,0)
#################Plot Survival TCGA###################
TCGA$cluster<- factor(TCGA$cluster, labels = c("C1", "C2", "C3"))
sfit1 <- survfit(Surv(days, vital_status) ~ cluster, TCGA)
plot(sfit1, mark.time=F, col=c(1,2,4), lty=1, lwd=2,
     xlab="Days from Sample",
     ylab="Survival")
legend("topleft", legend=levels(TCGA$cluster), text.col=c(1,2,4))
#################Plot Survival JP###################
JP$cluster<- factor(JP$cluster, labels = c("C1", "C2", "C3"))
sfit2 <- survfit(Surv(days, vital_status) ~ cluster,  JP)
plot(sfit2, mark.time=F, col=c(1,2,4), lty=1, lwd=2,
     xlab="Days from Sample",
     ylab="Survival")
legend("topleft", legend=levels(TCGA$cluster), text.col=c(1,2,4))


