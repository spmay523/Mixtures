library(methods)

EMmix <- setClass(
  "EMmix",
  slots = c(
    pnames = "matrix", x="data.frame", dat="data.frame", pi_init="numeric",
    med = "matrix", mu = "matrix", x90 = "matrix"
  )
)

setMethod("initialize",
          signature = c(.Object="EMmix"),
          function(.Object, pn, x, d, p_init, med, mu, x90){
            .Object@pnames = pn
            .Object@x = x
            .Object@dat = d
            .Object@pi_init = p_init
            .Object@med = med
            .Object@mu = mu
            .Object@x90 = x90
            return(.Object)
          }
)

setMethod("plot",
          signature = c(x="EMmix",y="missing"),
          function(x,y){
            last_iter <- nrow(x@med)
            pops <- x@pnames[,1]
            for(p in pops){
              xmin <- min(x@med[,p])
              xmax <- max(x@med[,p])
              last <- x@med[last_iter,p]
              plot.default(x=x@med[,p], main=p, xlab="Iteration Number", ylab = "Median Mixture Parameter", ylim=c(xmin,xmax),type="p")
              abline(h=last, lty=2, col="firebrick1")
            }
          }
)

setMethod("show",
          signature = "EMmix",
          function(object){
            last_iter <- nrow(object@med)
            cat(last_iter," iterations conducted\n\nFinal Median Mixture Parameters:\n")
            print(sort(object@med[last_iter,], decreasing = TRUE))
            cat("\n\nFinal Mean Mixture Parameters:\n")
            print(sort(object@mu[last_iter,], decreasing = TRUE))
          })
