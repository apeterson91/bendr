#' Estimate Clustered Inhomgogenous poisson process intensity functions from grouped data using a normal kernel
#'
#' @param formula bendr formula see the details
#' @param benvo  \code{\link[rBenvo]{Benvo}} object 
#' @param base_measure one of c("Normal","Beta") indicating which base measure to use.
#' @param L component truncation number
#' @param K intensity cluster truncation number
#' @param mu_0 mean hyperparameter for mu normal base measure. Default is 0; Normal(0,1).
#' @param kappa_0 variance hyperparameter for mu normal base measure. Default is 1; Normal(0,1).
#' @param nu_0 df hyperparameter for sigma inv chisq base measure. Default is 1; InvChisq(1,1);
#' @param sigma_0 scale hyperparameter for sigma inv chisq base measure. Default is 1; InvChisq(1,1);
#' @param a_alpha shape hyperparameter for alpha gamma prior
#' @param b_alpha scale hyperparameter for alpha gamma prior
#' @param a_rho shape hyperparameter for rho gamma prior
#' @param b_rho scale hyperparameter for rho gamma prior
#' @param iter_max total number of iterations for which to run sampler
#' @param burn_in number of iterations to discard as sampler approaches the stationary distribution
#' @param thin number of iterations to thin by
#' @param seed integer with which to initialize random number generator
#'
#' @export
bend <- function(formula,
				 benvo,
				 base_measure = "Normal",
				 mu_0 = 0, kappa_0 = 1,
				 nu_0 = 1, sigma_0 = 1,
				 L = 4 , K = 4,
				 a_alpha = 1, b_alpha = 1,
				 a_rho = 1, b_rho = 1,
				 iter_max = 5E3,
				 burn_in = floor(iter_max/2),
				 thin = 1,
				 seed = NULL) 
{

    call <- match.call(expand.dots=TRUE)
	dt <- groupvo(benvo,formula)
    ## basic checks
	stopifnot((iter_max > warm_up) && (warm_up > 0))
	stopifnot(L>0 && K >0)

    if(is.null(seed))
        seed <- 134143L

    d <- seq(from = 0, to = 1, by = 0.01) ## distance grid
    num_posterior_samples <- sum(seq(from=warm_up+1,to = iter_max,by=1) %% thin == 0 )
	stopifnot(num_posterior_samples >0 )


    fit <- nd_nhpp_fit(r= dt$r, n_j = dt$n, d = d,
					   L = L, K = K, J = dt$J,
					   mu_0 = mu_0, kappa_0 = kappa_0,
					   nu_0 = nu_0, sigma_0 = sigma_0,
					   a_alpha = a_alpha, b_alpha = b_alpha,
					   a_rho = a_rho, b_rho = b_rho,
					   iter_max = iter_max, warm_up = warm_up,
					   thin = thin, seed = seed, chain = 1,
					   num_posterior_samples = num_posterior_samples)


    d <- stats::pnorm(d)

	obj <- list(fit = fit,
				J = J,
				d = dt$R*d,
				call = call,
				K = K, L = L,
				chains = 1)

    out <- ndp(obj)
}

