setwd("/home/xnestea/Master/bin/Merged_MAD/")

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

install.packages("NMF")
library(NMF)

load("Exp_TCGA.Rdata")
model <- nmf(exp, 3, nrun = 200, .options = "v", .pbackend="mpi")
#consensusmap(estim.r)

save(model, file = "model_TCGA.Rdata")