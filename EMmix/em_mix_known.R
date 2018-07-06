########### Old code
###########  EM1 FUNCTION  #############
N=nrow(x) # number of markers

p<-exac[which(x.touse$AF>MAF_thresh & x.touse$AF<(1-MAF_thresh)),paste(c("AC_hom_", "AC_het_", "AC_homref_"), rep(pop_names, each=3), sep="")]
k=length(pop_names)

pi_out<-pi_out_median<-pi_new<-pi_start
iter=0
thresh_check=threshold+1

while( thresh_check > threshold){
  pi=pi_new
  gamma_tmp.a<-matrix(NA, nrow=nrow(x), ncol=k)
  for(p_k in 1:k){
    ##estimating the posterior probs from the data per k ancestry group per n markers
    x.p<-cbind(x, p[,c(((p_k-1)*3+1):((p_k-1)*3+3))])
    gamma_tmp.a[,p_k]<-apply(x.p,1, function(x.f){
      p_k_geno<-as.numeric(x.f[4:6])
      pi[p_k]*dmultinom(as.numeric(x.f[1:3]), prob = p_k_geno, log=T)
    })
  }
  
  gamma_tmp<-(apply(gamma_tmp.a,1,function(x){sum(x, na.rm=T)})-gamma_tmp.a)/apply(gamma_tmp.a,1,function(x){sum(x, na.rm=T)})
  
  pi_new<-apply(gamma_tmp,2,function(x){mean(x, na.rm=T)})
  pi_median<-apply(gamma_tmp,2,function(x){median(x, na.rm=T)})
  pi_90<-apply(gamma_tmp,2,function(x){quantile(x, probs=0.90,na.rm=T)})
  
  pi_out<-rbind(pi_out, pi_new)
  pi_out_median<-rbind(pi_out_median, pi_median)
  iter=iter+1
  thresh_check<-sum(abs(pi-pi_new))
  if(round(iter/10)==(iter/10)){write.table(gamma_tmp, paste(person, "/gamma_values/pops", pop_names[1], N_pop1, "_", pop_names[2], (N-N_pop1), "_start", paste(pi_start, collapse="_"), "_iter", iter, ".txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")}
  print(pi)
  print(nrow(x))
}

#' @param x N by K data.frame for N markers and K ancestries giving allele frequencies for each ancestry. Must have column named 'AF' for total allele frequency
#' @param MAF_thresh numeric vector of length 1 with value between 0 and 0.5 indicating when to exclude a marker due to rarity
#' @param pnames K by 4 character matrix of population names in position 1, followed by hom, het, and homref indicating the column names for each category
#' @param pi_init numeric vector of length K providing initial mixture parameters
#' @param threshold numeric vector of length 1 indicating when to terminate the estimation procedure
#' @param dat data.frame object for frequencies by each population
#' @param path gives a file path name for the directory where tables of outputs from the iterations should be written
#' @param Ntot total number of population
#' @param N_pop1 number of population 1
em_mix_known <- function(x, dat, MAF_thresh = 0.05, pnames, pi_init, threshold = 0.01, path = NA, Ntot, N_pop1){
  pi <- pi_new <- pi_median <- pi_init
  names(pi) <- pnames[,1]
  k <- nrow(pnames)
  N <- nrow(x)
  thresh_check <- threshold + 1
  cols_toSelect <- character()
  apply(pnames[,2:4], 1, function(pop){cols_toSelect <<- c(cols_toSelect, pop)} )
  p <- dat(which(x$AF > MAF_thresh & x$AF < (1-MAF_thresh)))
  iter = 0
  
  # expectation step
  while(thresh_check > threshold){
    pi = pi_new
    #new create gamma_tmp
    gamma_tmp.a <- cbind(
      apply(pnames, 1, function(pop){
        x.p <- cbind(x, p[,pop[2:4]])
        print(x.p)
        apply(x.p, 1, function(gen){
          pi[pop[1]]*dmultinom(as.numeric(gen[1:3]), prob = gen[4:6], log = T)
        })
      })
    )
    gamma_tmp<-(apply(gamma_tmp.a,1,function(x){sum(x, na.rm=T)})-gamma_tmp.a)/apply(gamma_tmp.a,1,function(x){sum(x, na.rm=T)})
    
    pi_new<-apply(gamma_tmp,2,function(x){mean(x, na.rm=T)})
    names(pi_new) <- names(pi)
    pi_median<-apply(gamma_tmp,2,function(x){median(x, na.rm=T)})
    pi_90<-apply(gamma_tmp,2,function(x){quantile(x, probs=0.90,na.rm=T)})
    
    pi_out<-rbind(pi_out, pi_new)
    pi_out_median<-rbind(pi_out_median, pi_median)
    iter=iter+1
    thresh_check<-sum(abs(pi-pi_new))
    if(!path){
      if(round(iter/10)==(iter/10)){
        write.table(gamma_tmp, paste(path,"/pops", names(pi)[1], N_pop1, "_", names(pi)[2], (N-N_pop1), "_start", paste(pi_start, collapse="_"), "_iter", iter, ".txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")
      }
    }
  }
  return(pi)
}