#' Estimates the mean number of observations in an Inhomogenous Poisson Process
#'
#' 
#' fits an exponential linear regression model with log link.
#' 
#' @template nuts
#'
#' @section Details: 
#' Fits an exponential linear regression model using a log link via
#' No-U-Turn Sampler (see reference). This model has a 
#'  \eqn{N(0,\sigma=3)} prior on all regression coefficients.
#' @param formula formula specifying the design matrix and outcome 
#' akin to \code{\link[stats]{glm}} for more information.
#' @param data data.frame from which to extract outcome and covariates
#' @param warm_up number of iterations in which to tune HMC step-size, these will be discarded
#' @param iter_max total number of samples for which to run sampler
#' @param seed random number generator intializing seed
#'
#' @export
nhpp_hmc <- function(formula,
					 data,
					 warm_up,
					 iter_max,
					 seed = NULL){

	stopifnot( (warm_up<iter_max) && (iter_max>0) && (warm_up >0) )
	if(is.null(seed))
		seed <- 23413

	call <- match.call(expand.dots = TRUE)
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula" ),
			 table = names(mf), nomatch = 0L)
	mf <- mf[c(1L, m)]
	mf$data <- data
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	Y <- stats::model.response(mf, type = "any")
	if (stats::is.empty.model(mt))
	stop("No intercept or predictors specified.", call. = FALSE)
	X <- stats::model.matrix(mt, mf)

	fit <- nhpp_gamma(warm_up,iter_max,X,Y,0.65,seed)

	return(hmc(fit,warm_up,iter_max,X))

}
