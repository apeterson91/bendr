#' Base Measures
#' 
#' Pros and Cons of different base measures.
#'
#' @details One of the "parameters" of the Dirichlet Process is the base measure
#' which describes the distribution around which the Dirichlet Process is centered.
#' \code{\link{bend}} allows for two possible base measures: (1) normal and (2) beta.
#' 
#' The normal base measure is faster and offers more efficient samples, but requires 
#' transforming the distances or radii to the real line via the nonlinear probit link function.
#' In contrast the beta base measure requires only a linear scaling to keep distances in (0,1)
#' and may be more appealing for it's lack of "edge effects" in the estimated densities.
#' 
#' @name measures
#' @seealso \link{normal_measure} \link{beta_measure}
#'
NULL


#' Normal Base Measure
#'
#' @param mu_0 mean hyperparameter for mu normal base measure. Default is 0; Normal(0,1).
#' @param kappa_0 variance hyperparameter for mu normal base measure. Default is 1; Normal(0,1).
#' @param nu_0 df hyperparameter for sigma inv chisq base measure. Default is 1; InvChisq(1,1);
#' @param sigma_0 scale hyperparameter for sigma inv chisq base measure. Default is 1; InvChisq(1,1);
#' @export
#' 
normal_measure <- function(mu_0 = 0,nu_0 = 1,sigma_0 = 1,kappa_0 = 1){

	structure(
			  list(measure = "normal",
					mu_0 = mu_0,
					nu_0 = nu_0,
					sigma_0 = sigma_0,
					kappa_0 = kappa_0),
					 class = "bm"
		  )
}


#' Beta base measure
#'
#' @param a_0 hyperparameter for mu beta base measure. Default is 1; uniform(0,1).
#' @param b_0 hyperparameter for mu beta base measure. Default is 1; uniform(0,1).
#' @param mu_sd proposal scale for mus 
#' @param tau_sd proposal scale for mus 
#' @export
#' 
beta_measure <- function(a_0 = 1,b_0 = 1,mu_sd = 1, tau_sd = 1){

	structure(
			list(measure = "beta",
				 a_0 = a_0,
				 b_0 = b_0,
				 mu_sd = mu_sd,
				 tau_sd = tau_sd), 
			class = "bm"
	)
}


#' Transform Distances 
#' 
#' @details Transforms distances according to which base measure is used.
#' Either a linear transformation - scaling by the ceiling of the maximum distance 
#' or a nonlinear transformation via the probit function.
#' @param bm base measure object
#' @param r vector of radii
#'
transform_distances <- function(bm,r)
	UseMethod("transform_distances")

transform_distances.bm <- function(bm,r){
	R <- ceiling(max(r))
	r_ <- r/R
	if(is_normal(bm))
		return(qnorm(r_))
	else
		return(return(r_))
}


is_normal <- function(bm)
	UseMethod("is_normal")

is_normal.bm <- function(bm){
	if(bm$measure=="normal")
		return(TRUE)
	else
		return(FALSE)
}
