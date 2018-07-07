#' @param x N by K data.frame for N markers and K ancestries giving allele frequencies for each ancestry. Must have column named 'AF' for total allele frequency
#' @param MAF_thresh numeric vector of length 1 with value between 0 and 0.5 indicating when to exclude a marker due to rarity
#' @param pnames matrix of population names in column 1 that contains K rows, followed by the column names pertaining to the allele frequency(ies) for that population as they appear in dat.
#' @param pi_init numeric vector of length K providing initial mixture parameters
#' @param threshold numeric vector of length 1 indicating when to terminate the estimation procedure
#' @param dat data.frame object for frequencies by each population
#' @param path gives a file path name for the directory where tables of outputs from the iterations should be written
#' @param Ntot total number of population
#' @param N_pop1 number of population 1
#' @return Numeric vector of mixtures for each population included as a row in pnames
#' @export
em_mix_known <- function(x, dat, pnames, pi_init, threshold = 0.01, MAF_thresh = 0.05, path = NA){
  pi_out <- pi_out_median <- pi_out_90 <- pi_new <- pi_median <- pi_90 <- pi_init
  names(pi_median) <- names(pi_new) <- names(pi_90) <- pnames[,1]
  k <- nrow(pnames)
  N <- nrow(x)
  thresh_check <- threshold + 1
  cols_toSelect <- character(); apply(pnames[,2:ncol(pnames)], 1, function(pop){cols_toSelect <<- c(cols_toSelect, pop)} )
  p <- dat # p <- dat(which(x$AF > MAF_thresh & x$AF < (1-MAF_thresh)))
  iter = 0
  
  # expectation step
  while(thresh_check > threshold){
    pi = pi_new
    #new create gamma_tmp
    gamma_tmp.a <- cbind(
      apply(pnames, 1, function(pop){
        x.p <- cbind(x, p[,pop[2:ncol(pnames)]])
        apply(x.p, 1, function(gen){
          pi[pop[1]]*dmultinom(as.numeric(gen[1:(ncol(pnames)-1)]), prob = gen[ncol(pnames):(2 * (ncol(pnames)-1))], log = T)
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
    pi_out_90<-rbind(pi_out_90, pi_90)
    iter=iter+1
    thresh_check<-sum(abs(pi-pi_new))
    # the write to file is currently not implemented; this will likely need to be deleted
    # if(!is.na(path)){
    #   if(round(iter/10)==(iter/10)){
    #     write.table(gamma_tmp, paste(path,"/pops", names(pi)[1], N_pop1, "_", names(pi)[2], (N-N_pop1), "_start", paste(pi_start, collapse="_"), "_iter", iter, ".txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")
    #   }
    # }
  }
  iter_results <- cbind(pi_out, pi_out_median, pi_out_90)
  dimnames(iter_results) <- list(c(),c(
    c(paste(pnames[,1], "_mean", sep = "")),
    c(paste(pnames[,1], "_median", sep = "")),
    c(paste(pnames[,1], "_90", sep = ""))
  ))
  print(pi_median)
  return(iter_results)
}