#' Print method for ndp objects
#' 
#' The \code{print} method for ndp objects displays a compact summary of the
#' fitted model. See the \strong{Details} section below for descriptions of the
#' different components of the printed output. For additional summary statistics
#' and diagnostics use the \code{\link[=summary.ndp]{summary}} method.
#' 
#' @export
#' @param x a ndp object
#' @param digits Number of digits to use for formatting numbers.
#' @param ... Ignored.
#' @return Returns \code{x}, invisibly. 
#' @details 
#' \subsection{Point Estimates}{
#' Point estimates are medians computed from simulations.
#' }
#' \subsection{Uncertainty estimates (SD)}{
#' The standard deviations reported (labeled \code{MAD} in the print output)
#' are computed from the same set of draws described above and are proportional
#' to the median absolute deviation (\code{\link[stats]{mad}}) from the median.
#' Compared to the raw posterior standard deviation, the MAD_SD will be
#' more robust for long-tailed distributions. These are the same as the values
#' }
#' 
#' @seealso \code{\link{summary.ndp}}
#' @importFrom coda as.matrix
#' 
print.ndp <- function(x, digits = 2, ...) {
  cat("\n observations:", x$n)
  cat("\n groups:", x$J)

  cat("\n------\n")
  cat("Cluster Statistics")
  cat("\n------\n")

  pi_stats <- rbind(Median = apply(as.matrix(x$pi),2,stats::median),
                    MAD = apply(as.matrix(x$pi),2,stats::mad))

  .printfr(pi_stats,digits)

  if(all(as.matrix(x$alpha) == 0 ))
      cat("alpha collapsed at 0 \n")
  else{
      .printfr(cbind(Median = apply(as.matrix(x$alpha),2,stats::median),
                MAD = apply(as.matrix(x$alpha),2,mad)), digits)
      cat("\n")
      }
  if(all(as.matrix(x$rho) == 0 ))
      cat("rho collapsed at 0 " )
  else{
      .printfr(cbind(Median = apply(as.matrix(x$rho),2,stats::median),
                MAD = apply(as.matrix(x$rho),2,mad)),digits)
      cat("\n")
  }    

  cat(" \n---\n ") 

  mat <- cbind(min = min(sapply(x$num_clusters,min)),
            median = stats::median(sapply(x$num_clusters,stats::median)),
            max = max(sapply(x$num_clusters,max)))
  rownames(mat) <- "# of Clusters"

  .printfr(mat,digits)
  cat("\n------\n")

  w_stats <- rbind(Median = apply(as.matrix(x$w),2,stats::median),
                    MAD = apply(as.matrix(x$w),2,mad))
  mu_stats <- rbind(Median = apply(as.matrix(x$mu),2,stats::median),
                    MAD = apply(as.matrix(x$mu),2,mad))
  tau_stats <- rbind(Median = apply(as.matrix(x$tau),2,stats::median),
                     MAD = apply(as.matrix(x$tau),2,mad))

  tau_switch <- all(dim(tau_stats) == dim(mu_stats))

  K <- max(x$if_df$Intensity_Function)
  L <- ncol(w_stats)/K

  for(k in 1:K){
      start <- (1 + L*(k-1))
      end <- (L + L*(k-1))
      if(pi_stats["Median",k]>.01){
          .printfr(w_stats[,start:end],digits)
          .printfr(mu_stats[,start:end],digits)
          if(tau_switch)
              .printfr(tau_stats[,start:end],digits)
      }
      cat("\n")
  }
  if(!tau_switch)
      .printfr(tau_stats,digits)

  
  cat("\n------\n")
  cat("* For help interpreting the printed output see ?print.ndp\n")
  
  invisible(x)
}


#' Summary method for ndp objects
#' 
#' Summaries of parameter estimates and MCMC convergence diagnostics 
#' (Monte Carlo error, effective sample size, Rhat).
#' 
#' @export
#' 
#' 
#' @param object  an ndp object
#' @param ... Currently ignored.
#' @param digits Number of digits to use for formatting numbers when printing. 
#'   When calling \code{summary}, the value of digits is stored as the 
#'   \code{"print.digits"} attribute of the returned object.
#'   
#' @return The \code{summary} method returns an object of class 
#'   \code{"summary.ndp"}.
#' 
#' 
#' @importFrom coda effectiveSize 
#' @importFrom stats sd quantile
summary.ndp <- function(object,digits = 1, ...) {


    summ <- function(x,y){ 
        out <- t(apply(as.matrix(x[[y]]),2,function(a) c("Mean"=mean(a),"SD" = sd(a), quantile(a),"Geweke" = unname(coda::geweke.diag(a)$z) )))
        out <- cbind(out,"ESS"=effectiveSize(x[[y]]))
        if(length(x[[y]])>1)
            out <- cbind(out,"Rhat" = coda::gelman.diag(x[[y]],multivariate = FALSE,autoburnin=FALSE)$psrf[,1] )
        return(out)
        }
    # alpha and rho
    out <- rbind(alpha = summ(object,"alpha"),
                 rho = summ(object,"rho"))
    out <- rbind(out,summ(object,"pi"),summ(object,"w"),summ(object,"mu"),summ(object,"tau"))


  structure(
    out,
    nobs = object$n,
    groups = object$J,
    posterior_sample_size = nrow(as.matrix(object$alpha)),
    call = object$call,
    print.digits = digits,
    class = "summary.ndp"
  )

}


#' @rdname summary.ndp
#' @export
#' @method print summary.ndp
#'
#' @param x An object of class \code{"summary.ndp"}.
print.summary.ndp <- function(x, digits = max(1, attr(x, "print.digits")), 
                                  ...) {
    atts <- attributes(x)
    cat("\nModel Info:\n")
  cat("\n sample:      ", atts$posterior_sample_size, "(posterior sample size)")
  cat("\n observations:", atts$nobs)
  cat("\n groups:", atts$groups)

  cat("\n\nEstimates:\n")
    .printfr(x,digits)

    invisible(x)
}

#' @export
print.hmc <- function(x, digits = 2, ...) {
  cat("\n groups:", x$J)

  cat("\n------\n")
  cat("Regression Statistics")
  cat("\n------\n")


  if(dim(x$beta)[2]==1){
      meds <- summary(x$beta)$quantile["50%"]
      sd <- summary(x$beta)$statistics["SD"]
    }
  if(dim(x$beta)[2]>1){
      meds <- summary(x$beta)$quantile[,"50%"]
      sd <- summary(x$beta)$statistics[,"SD"]
    }

  mat <- cbind(Median = meds, SD = sd)
  rownames(mat) <- colnames(x$X)

    .printfr(mat, digits, ...)

  cat("--- \n") 

  ppd <- coda::as.mcmc(apply(x$beta,1,function(beta) mean(exp((x$X %*% beta) )) ))
  mat <- cbind(Median = summary(ppd)$quantile["50%"],
               SD = summary(ppd)$statistics["SD"])
  rownames(mat) <- "mean_PPD"
  cat("\nSample avg. posterior predictive distribution of y:\n")
  .printfr(mat,digits,...)
  
  cat("\n------\n")
  
  invisible(x)
}


# internal ----------------------------------------------------------------
.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
}
