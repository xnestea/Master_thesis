setwd("/home/xnestea/Master/bin")

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

install.packages("NMF")
library(NMF)

load("TPM/prep_exp_TCGA.Rdata")
model <- nmf(exp, 2, nrun = 200, .options = "v", .pbackend="mpi")


save.image(file = "model_TCGA.Rdata")