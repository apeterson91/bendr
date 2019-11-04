#' Estimate the nonhomgogenous poisson process intensity function from grouped data using a beta base measure
#' 
#' @param r vector of distances  associated with differing groups
#' @param n_j matrix of integers denoting the start and length of each school's associated BEF distances
#' @param mu_sd proposal scale for mus 
#' @param tau_sd proposal scale for mus 
#' @param L component truncation number
#' @param K intensity cluster truncation number
#' @param a_0 hyperparameter for mu beta base measure. Default is 1; uniform(0,1).
#' @param b_0 hyperparameter for mu beta base measure. Default is 1; uniform(0,1).
#' @param a_alpha hyperparameter for alpha gamma prior
#' @param b_alpha hyperparameter for alpha gamma prior
#' @param a_rho hyperparameter for rho gamma prior
#' @param b_rho hyperparameter for rho gamma prior
#' @param iter_max total number of iterations for which to run sampler
#' @param warm_up number of iterations for which to burn-in or "warm-up" sampler
#' @param thin number of iterations to thin by
#' @param seed integer with which to initialize random number generator
#'
#' @export
beta_nd_nhpp <- function(r, n_j,
                    mu_sd, tau_sd, 
                    L = 4 , K = 4, 
                    a_0 = 1 , b_0 = 1, 
                    a_alpha = 1, b_alpha = 1,
                    a_rho = 1, b_rho = 1, 
                    iter_max, warm_up, 
                    thin = 1,
                    multiple_taus = FALSE,
                    kernel = "Normal",
                    seed = NULL) {

    call <- match.call(expand.dots=TRUE)
    J <-  nrow(n_j)
    ## basic checks
    if(any(r<=0))
        stop("all r must be positive numbers")
    if(any(r>=1)){
        R <- ceiling(max(r)) ## maximum radius
        r_ <- r / R ## scaled event radii 
    }else{
        r_ <-  r
        R <- 1
    }
    
    d <- seq(from = 0.01, to = 0.99, by = 0.01) ## distance grid

    if(is.null(seed))
        seed <- 1L

	if(multiple_taus)
		fit <- list(beta_nd_nhpp_fit_multiple_taus(r = r_,n_j = n_j,d = d,
						 mu_sd = mu_sd, tau_sd = tau_sd,
						 L =  L,K = K, J = J,
						 a_0 = a_0, b_0 = b_0,
						 a_alpha = a_alpha ,b_alpha = b_alpha,
						 a_rho = a_rho, b_rho = b_rho,
						 iter_max = iter_max, warm_up = warm_up,
						 thin = thin,seed = seed, chain = 1L))
	else
		fit <- list(beta_nd_nhpp_fit( r = r_,n_j = n_j,d = d,
											  mu_sd = mu_sd, tau_sd = tau_sd,
											  L =  L,K = K, J = J,
											  a_0 = a_0, b_0 = b_0,
											  a_alpha = a_alpha ,b_alpha = b_alpha,
											  a_rho = a_rho, b_rho = b_rho,
											  iter_max = iter_max, warm_up = warm_up,
											  thin = thin,seed = seed, chain = 1L))
            

    out <- bndp(c(list(K = K, L = L, d = R*d, n = sum(n_j[,2]), call = call),fit),1,J)
}

