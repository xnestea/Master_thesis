setwd("/home/xnestea/Master/bin")

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

install.packages("NMF")
library(NMF)

load("../bin/prep_exp_TCGA.Rdata")
estim.r <- nmf(exp, 2:8, nrun = 50, .options = "v", .pbackend="mpi")
plot_coef <- plot(estim.r)
#consensusmap(estim.r)

save.image(file = "k_estimation_TCGA.Rdata")
