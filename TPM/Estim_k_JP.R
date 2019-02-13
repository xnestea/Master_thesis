setwd("/home/xnestea/Master/bin")

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

install.packages("NMF")
library(NMF)

load("../bin/prep_exp_JP.Rdata")
estim.r <- nmf(exp, 2:8, nrun = 50, .options = "v", .pbackend="mpi")
plot_coef <- plot(estim.r)

save.image(file = "k_estimation_JP.Rdata")