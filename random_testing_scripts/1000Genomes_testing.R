setwd("/home/projects/mixtures/Sam/EMmix_repo/Mixtures/EMmix/")
source("em_mix_known.R")
source("em_object.R")

dat <- read.table("/home/projects/mixtures/yinfei/temp4", header = TRUE)
popNames <- c("AMR", "AFR", "EUR", "EAS", "SAS")
l <- nrow(dat)
hardyWeinNames <- paste(rep(popNames, each = 3), c("_hom","_het","_homref"),sep="")
dat_hardyWein <- as.data.frame(matrix(nrow=l, ncol=length(hardyWeinNames), dimnames = list(dat$POS,hardyWeinNames)))
pnames <- matrix(nrow = length(popNames), ncol=4, dimnames = list(popNames, c()))
for(pop in popNames){
  column <- dat[,paste(pop, "_AF", sep = "")]
  # assume hardy-weinberg equilibrium
  hom <- column ^ 2
  het <- 2 * column * (1-column)
  homref <- (1-column)^2
  cols <- paste(rep(pop, each = 3), c("_hom","_het","_homref"),sep="")
  dat_hardyWein[,cols] <- cbind(hom, het, homref)
  pnames[pop,] <- c(pop, cols)
}

num_pop <- 1000
real_pi <- c(0.3,0,0.2,0,0.5)
pop_number <- num_pop  * real_pi
names(pop_number) <- popNames
names(real_pi) <- popNames
pop_sim<-as.data.frame(matrix(nrow = l, ncol = length(hardyWeinNames), dimnames = list(dat$POS,hardyWeinNames)))
for(pop in popNames){
  popDat <- dat_hardyWein[,paste(rep(pop, each = 3), c("_hom","_het","_homref"),sep="")]
  pop_sim <<- cbind(pop_sim, t(apply(popDat, 1, function(pd){rmultinom(pop_number, size = 1, prob = pd)})))
}

# next need to get the overall frequency from the simulated columns
sim_x <- data.frame(hom = apply(pop_sim[,which(1:(3 * length(popNames)) %% 3 == 1)], 1, sum),
                    het = apply(pop_sim[,which(1:(3 * length(popNames)) %% 3 == 2)], 1, sum),
                    homref = apply(pop_sim[,which(1:(3 * length(popNames)) %% 3 == 0)], 1, sum)
                   )
# create pnames for referencing dat_hardyWein
mixture <- em_mix_known(x = sim_x, dat = dat_hardyWein, pnames = pnames, pi_init = c(0.2,0.2,0.2,0.2,0.2), Ntot = num_pop)