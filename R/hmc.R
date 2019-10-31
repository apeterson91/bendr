#' Create a HMC model object
#' 
#' @param object model fit object returned from nhpp_gamma Rcpp function
#' @param warm_up same as nhpp_hmc
#' @param iter_max same as nhpp_hmc
#' @param X the design matrix used in the regression
#' @return an hmc model object
hmc <- function(object,warm_up,iter_max,X){

	ics <- intersect((warm_up+1):iter_max,which(object$acceptance==1))
	
	beta <- coda::as.mcmc(object$beta_samples[ics,,drop=F])
	
	return(structure(list(beta = beta,
						  treedepth = object$treedepth,
						  epsilons = object$epsilons,
						  X = X),class = c("hmc")))

}
