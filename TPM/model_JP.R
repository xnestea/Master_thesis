setwd("/home/xnestea/Master/bin")

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

install.packages("NMF")
library(NMF)
load("../bin/prep_exp_JP.Rdata")
### in case there are full null rows remove them
exp <- exp[rowSums(exp)!= 0, ]
model <- nmf(exp, 2, nrun = 200, .options = "v", .pbackend="mpi")

save.image(file = "model_JP.Rdata")