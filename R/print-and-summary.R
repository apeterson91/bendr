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
#'
print.ndp <- function(x, digits = 2, ...) {
  cat("\n Base Measure:", x$base_measure$measure)
  cat("\n formula:", formula_string(formula(x$call)))
  cat("\n observations:", x$n)
  cat("\n groups:", x$J)

  cat("\n------\n")
  cat("Cluster Statistics")
  cat("\n------\n")

  pi_stats <- rbind(Median = apply(x$p,2,stats::median),
                    MAD = apply(x$p,2,stats::mad))

  .printfr(pi_stats,digits)


  cat(" \n---\n ")
  nc <- num_clusters(x$cluster_assignment)

  mat <- cbind(min = min(nc),
			   median = median(nc),
			   max = max(nc))
  rownames(mat) <- "# of Clusters"

  .printfr(mat,digits)
  cat("\n------\n")


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
#' @importFrom stats sd quantile
summary.ndp <- function(object,digits = 1, ...) {

	parmat <- as.matrix(object)
	mean <- colMeans(parmat)
	sd <- apply(parmat,2,sd)
	qs <- t(apply(parmat,2,function(x) quantile(x,c(0.1,.25,.5,.75,.9),na.rm=T)))
	n_eff <- apply(parmat,2,rstan::ess_tail)
	Rhat <- apply(parmat,2,rstan::Rhat)
	out <- cbind(mean,sd,qs,n_eff,Rhat)


  structure(
    out,
    nobs = object$n,
	bm = object$base_measure,
    groups = object$J,
    posterior_sample_size = nrow((object$alpha)),
	formula = formula(object$call),
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
print.summary.ndp <- function(x,
							  digits = max(1, attr(x, "print.digits")),
							  ...) {
  atts <- attributes(x)
  cat("\nModel Info:\n")
  cat("\n Base Measure:", atts$bm$measure)
  cat("\n sample:      ", atts$posterior_sample_size, "(posterior sample size)")
  cat("\n observations:", atts$nobs)
  cat("\n groups:", atts$groups)

	hat <- "Rhat"
	str_diag <- "MCMC diagnostics"
	str1 <- "and Rhat is the potential scale reduction factor on split chains"
	str2 <- " (at convergence Rhat=1).\n"

	sel <- setdiff(colnames(x),c("n_eff","Rhat"))
	xtmp <- x[,sel]
	cat("\n\n")

    # print table of parameter stats
    .printfr(xtmp, digits)

	cat("\n", str_diag, "\n", sep = '')
	mcse_hat <- format(round(x[, c(hat), drop = FALSE], digits),
					  nsmall = digits)
	n_eff <- format(x[, "n_eff", drop = FALSE], drop0trailing = TRUE)
	print(cbind(mcse_hat, n_eff), quote = FALSE)
	cat("\n n_eff is a crude measure of effective sample size, ",
	  str1,
	  str2, sep = '')

    invisible(x)
}


# internal ----------------------------------------------------------------
.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
}

num_clusters <- function(cluster_mat){

	  apply(cluster_mat,1,function(x) length(unique(x)))

}

formula_string <- function(formula, break_and_indent = TRUE) {
  coll <- if (break_and_indent) "--MARK--" else " "
  char <- gsub("\\s+", " ", paste(deparse(formula), collapse = coll))
  if (!break_and_indent)
    return(char)
  gsub("--MARK--", "\n\t  ", char, fixed = TRUE)
}

.median_and_madsd <- function(x) {
  cbind(Median = apply(x, 2, median), MAD_SD = apply(x, 2, mad))
}
