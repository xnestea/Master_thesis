setwd("/home/xnestea/Master/bin/Merged_MAD")

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

install.packages("NMF")
library(NMF)

load("Exp_TCGA.Rdata")
estim.r <- nmf(exp, 2:8, nrun = 50, .options = "v", .pbackend="mpi")
#consensusmap(estim.r)

save(estim.r, file = "k_estimation_TCGA.Rdata")
