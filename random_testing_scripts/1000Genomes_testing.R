setwd("/home/projects/mixtures/Sam/EMmix_repo/Mixtures/EMmix/")
source("em_mix_known.R")
source("em_object.R")

dat <- read.table("/home/projects/mixtures/yinfei/temp4")
popNames <- c("AMR", "AFR", "EUR", "EAS", "SAS")
l <- nrow(dat)
hardyWeinNames <- paste(rep(popNames, each = 3), c("_hom","_het","_homref"),sep="")
dat_hardyWein <- as.data.frame(matrix(nrow=l, ncol=length(hardyWeinNames), dimnames = list(dat$POS,hardyWeinNames)))
for(pop in popNames){
  column <- dat[,paste(pop, "_AF", sep = "")]
  # assume hardy-weinberg equilibrium
  hom <- column ^ 2
  het <- 2 * column * (1-column)
  homref <- (1-column)^2
  cols <- paste(rep(pop, each = 3), c("_hom","_het","_homref"),sep="")
  dat_hardyWein[,cols] <- cbind(hom, het, homref)
}

num_pop <- 1000
real_pi <- c(0.3,0,0.2,0,0.5)
pop_number <- num_pop  * real_pi
names(pop_number) <- popNames
names(real_pi) <- popNames
pop_sim<-as.data.frame(matrix(nrow = l))
for(pop in popNames){
  popDat <- dat_hardyWein[,paste(rep(pop, each = 3), c("_hom","_het","_homref"),sep="")]
  pop_sim <<- cbind(pop_sim, t(apply(popDat, 1, function(pd){rmultinom(pop_number, size = 1, prob = pd)})))
}
# next need to get the overall frequency from the simulated columns