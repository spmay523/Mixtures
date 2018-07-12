library(methods)

EMmix <- setClass(
  "EMmix",
  slots = c(
    pnames = "matrix", x="data.frame", dat="data.frame", pi_init="numeric",
    em_iter="matrix"
  )
)

setMethod("initialize",
          signature = c(.Object="EMmix"),
          function(.Object, pn, x, d, p_init, results){
            .Object@pnames = pn
            .Object@x = x
            .Object@dat = d
            .Object@pi_init = p_init
            .Object@em_iter = results
            return(.Object)
          }
)

setMethod("plot",
          signature = c(x="EMmix",y="missing"),
          function(x,y){
            last_iter <- nrow(x@em_iter)
            pops <- x@pnames[,1]
            for(p in pops){
              r <- x@em_iter[,paste(p,"_median",sep="")]
              xmin <- min(r)
              xmax <- max(r)
              last <- x@em_iter[last_iter,paste(p,"_median",sep="")]
              plot.default(x=r, main=p, xlab="Iteration Number", ylab = "Median Mixture Parameter", ylim=c(xmin,xmax),type="p")
              abline(h=last, lty=2, col="firebrick1")
            }
          }
)
