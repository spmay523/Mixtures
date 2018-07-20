# Running em_mix_known sample
setwd("/home/projects/mixtures/Sam/EMmix_repo")
source("Mixtures/EMmix/em_mix_known.R") # load EMmix object definition
source("Mixtures/EMmix/em_mix_known.R") # load em_mix_known function
Ntot = 5000
N = 1500 # this just limits the number of rows being pulled in; if you want to use the full data set then, remove the nrows parameter from the read.table command below
exac <- read.table("/home/projects/mixtures/ExAC/all.common_added_cols_v2.txt", sep="\t", header=T, as.is=T, fill=T, nrows=N)
# set up pnames key for populations of interest
pops <- c("AFR", "NFE", "SAS", "EAS")
pnames <- matrix(nrow = length(pops), ncol = 2)
for(i in 1:length(pops)){
  pnames[i,] <- c(pops[i], paste("AF_",pops[i],sep=""))
}
# simulate some data for allele frequencies
real_mix <- c(0.1, 0.3, 0.2, 0.3)
x <- t(rbind(apply(exac[,pnames[,2]], 1, function(x){
  sapply(x, FUN = function(i){rbinom(n = 1, size = Ntot, prob=i)})
})))

pi_init <- c(0.1, 0.2, 0.3, 0.4)

results <- em_mix_known(x = x, dat = exac, pnames = pnames, pi_init = pi_init)

write.table(results, file = "../em_mix_known_exacResults.txt", quote = FALSE, sep = "\t")