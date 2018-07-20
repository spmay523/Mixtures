#### test data for em_mix_known; easy and quick to run on a local device
setwd("C:/Users/spmay/github/EMmix_repo/Mixtures/EMmix")
source("em_mix_known.R")
source("em_object.R")
n <- 1000 #set number of markers to simulate
p <- rmultinom(n, 5000, c(0.1, 0.1, 0.1, 0.2, 0.2, 0.3))
p <- t(p)
pnames <- rbind(c("AFR",paste(c("AC_hom_", "AC_het_", "AC_homref_"), rep("AFR", each=3), sep = "")),
             c("NFE",paste(c("AC_hom_", "AC_het_", "AC_homref_"), rep("NFE", each=3), sep = ""))
            )
x <- t(rmultinom(n, 5000, c(0.2, 0.3, 0.5)))
names_p <- character()
apply(pnames[,2:4], 1, function(x){names_p <<- c(names_p,x)})
colnames(p) <- names_p

results <- em_mix_known(x = x, dat=p, pnames = pnames, pi_init = c(0.1,0.9))
plot(results)

