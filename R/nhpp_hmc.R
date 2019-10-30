#' Estimates the mean number of observations in an Inhomogenous Poisson Process
#'
#' @param formula 
#' @param data data.frame from which to extract outcome and covariates
#' @param warm_up number of iterations in which to tune HMC step-size, these will be discarded
#' @param iter_max total number of samples for which to run sampler
#' @param seed
#'
#' @export
nhpp_hmc <- function(formula,
					 data,
					 warm_up,
					 iter_max,
					 seed = NULL){

	
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
  Y <- model.response(mf, type = "any")
  if (is.empty.model(mt))
    stop("No intercept or predictors specified.", call. = FALSE)
  X <- model.matrix(mt, mf)

  fit <- nhpp_gamma(warm_up,iter_max,X,Y,0.65,seed)

  return(fit)

}
