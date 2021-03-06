#######################################
##  plots for EM-1 output            ##
##  will need to update code         ##
##  to make more flexible to         ##
##  different simulation parameters  ##
#######################################

##this code only works for two ancestries - will need to update for additional ancestries ##

person="Kendra" ## make a directory with your name - make sure you use your name when running the code so you are not overwriting others work

setwd(paste("/home/projects/mixtures", person, sep="/"))

int.list<-c(0.1, 0.5, 0.9) ##initial proportion values
sim.list<-c(0.1, 0.5, 0.9)  ## simulated proportion values for one ancestry
N=10000  ## simulated sample size

pops<-c('AFR', 'NFE') ##two populations simulated

##initializing plot##
png("AFR_NFE_only_100517_meanANDmedian.png", width=1500, height=500)
par(mai=c(0.5, 0.5, 0.1, 0.1), mfcol=c(3,3))

##looping through plotting parameters##
for(int in int.list){
for(sim in sim.list){
temp<-read.table(paste0("pops", pops[1], "_", sim*N, '_', pops[2], (N-sim*N), "_start", int, "_", (1-int), ".txt"), header=F, as.is=T)

plot(0:(nrow(temp)-1),temp[,1], ylim=c(0,1), pch=20, cex=2, cex.axis=2, cex.lab=2, xlab="", ylab="")
points(0,temp[1,1], cex=2, pch=8)
lines(0:(nrow(temp)-1),temp[,1], lty=1, lwd=2)
abline(h=sim, lty=1, cex=2)

points(0:(nrow(temp)-1),temp[,3], pch=18, cex=2, col="blue")
lines(0:(nrow(temp)-1),temp[,3], lty=3, lwd=2, col="blue")
}##end sim loop
} ##end int loop
dev.off()
