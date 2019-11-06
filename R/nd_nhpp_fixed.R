#' Estimate the nonhomgogenous poisson process intensity function from grouped data using a normal kernel
#'
#' @param r vector of distances  associated with differing groups
#' @param n_j matrix of integers denoting the start and length of each school's associated BEF distances
#' @param L component truncation number
#' @param K intensity cluster truncation number
#' @param mu_0 mean hyperparameter for mu normal base measure. Default is 0; Normal(0,1).
#' @param kappa_0 variance hyperparameter for mu normal base measure. Default is 1; Normal(0,1).
#' @param nu_0 df hyperparameter for sigma inv chisq base measure. Default is 1; InvChisq(1,1);
#' @param sigma_0 scale hyperparameter for sigma inv chisq base measure. Default is 1; InvChisq(1,1);
#' @param alpha outer DP concentration parameter
#' @param rho inner DP concentration parameter
#' @param iter_max total number of iterations for which to run sampler
#' @param warm_up number of iterations for which to burn-in or "warm-up" sampler
#' @param thin number of iterations to thin by
#' @param seed integer with which to initialize random number generator
#'
#' @export
nd_nhpp_fixed <- function(r, n_j,
						  mu_0 = 0, kappa_0 = 1,
						  nu_0 = 1, sigma_0 = 1,
						  alpha = 1, rho = 1,
						  L = 4 , K = 4,
						  iter_max, warm_up,
						  thin = 1,
						  seed = NULL) {

    call <- match.call(expand.dots=TRUE)
    J <-  nrow(n_j)
    ## basic checks
    if(iter_max <= warm_up)
      stop("warm_up must be < iter_max",.call = FALSE)
    if(any(r<=0))
        stop("all r must be positive numbers", .call = FALSE)
    if(any(r>=1)){
        R <- ceiling(max(r)) ## maximum radius
        r_ <- r / R ## scaled event radii
    }else{
        r_ <-  r
        R <- 1
    }



    if(is.null(seed))
        seed <- 1L

	r_  <-  stats::qnorm(r_)

    d <- seq(from = floor(min(r_)), to = ceiling(max(r_)), by = 0.01) ## distance grid
    num_posterior_samples <- sum(seq(from=warm_up+1,to = iter_max,by=1) %% thin == 0 )
    fit <- list(nd_nhpp_fixed_fit(r= r_, n_j = n_j, d = d,
								  L = L, K = K, J = J,
								  mu_0 = mu_0, kappa_0 = kappa_0,
								  nu_0 = nu_0, sigma_0 = sigma_0,
								  alpha = alpha, rho = rho,
								  iter_max = iter_max, warm_up = warm_up,
								  thin = thin, seed = seed, chain = 1,
								  num_posterior_samples = num_posterior_samples))
    d <- stats::pnorm(d)

    out <- ndp_fixed(c(list(K = K, L = L, d = R*d,
                      n = sum(n_j[,2]), call = call),fit),1,alpha,rho,J)
}

