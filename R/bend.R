#' Estimate the distribution of (B)uilt (E)nvironment amenities via the (N)ested (D)irichlet Process
#'
#' @param formula bendr formula see the details
#' @param benvo  \href{https://apeterson91.github.io/rbenvo}{\code{Benvo}} object
#' @param base_measure "bm" object from one of the base measure functions. See \link{measures}.
#' @param L component truncation number
#' @param K intensity cluster truncation number
#' @param a_alpha shape hyperparameter for alpha gamma prior. Fixed alpha if fix_concentration = TRUE.
#' @param b_alpha scale hyperparameter for alpha gamma prior
#' @param a_rho shape hyperparameter for rho gamma prior. Fixed rho if fix_concentration = TRUE.
#' @param b_rho scale hyperparameter for rho gamma prior
#' @param iter_max total number of iterations for which to run sampler
#' @param burn_in number of iterations to discard as sampler approaches the stationary distribution
#' @param thin number of iterations to thin by
#' @param fix_concentration boolean value indicating whether or not concentration parameters should be kept fixed. Currently only implimented for alpha base measure
#' @param seed integer with which to initialize random number generator
#'
#' @details Fits the Nested Dirichlet Process to data as specified in the formula and hyperparameter settings. Note that a 'bendr formula' simply has the grouping id on the left side of the tilde and the built environment feature name on the right hand side. 
#' No intercepts, or other terms are accepted or needed. 
#' @seealso The introductory \href{https://apeterson91.github.io/bendr/articles/Introduction.html}{vignette} for an example on how this function is used.
#'
#' @importFrom stats median quantile formula mad pnorm qnorm
#' @importFrom utils tail
#'
#' @export
bend <- function(formula,
				 benvo,
				 base_measure = normal_measure(),
				 L = 4 , K = 4,
				 a_alpha = 1, b_alpha = 1,
				 a_rho = 1, b_rho = 1,
				 iter_max = 5E3,
				 burn_in = floor(iter_max/2),
				 thin = 1,
				 fix_concentration = FALSE,
				 seed = NULL)
{

    call <- match.call(expand.dots=TRUE)
	dt <- groupvo(benvo,formula)
    ## basic checks
	stopifnot((iter_max > burn_in) && (burn_in > 0))
	stopifnot(L>0 && K >0)
	stopifnot(base_measure$measure %in% c("normal","beta"))
	if(base_measure$measure == "beta")
		stop("beta base measure currently not working")

    if(is.null(seed))
        seed <- 134143L


	r <- transform_distances(base_measure,dt$r)
    d <- seq(from = floor(min(r)), to = ceiling(max(r)), length.out = 300) ## distance grid
    num_posterior_samples <- sum(seq(from=burn_in+1,to = iter_max,by=1) %% thin == 0 )
	stopifnot(num_posterior_samples >0 )


	if(base_measure$measure == "normal"){
		if(!fix_concentration){
			fit <- nd_nhpp_fit(r= r,
							   n_j = dt$n,
							   d = d,
							   L = L, K = K, J = dt$J,
							   mu_0 = base_measure$mu_0, 
							   kappa_0 = base_measure$kappa_0,
							   nu_0 = base_measure$nu_0, 
							   sigma_0 = base_measure$sigma_0,
							   a_alpha = a_alpha, b_alpha = b_alpha,
							   a_rho = a_rho, b_rho = b_rho,
							   iter_max = iter_max, warm_up = burn_in,
							   thin = thin, seed = seed, chain = 1,
							   num_posterior_samples = num_posterior_samples)
		}else{
			fit <- nd_nhpp_fixed_fit(r= r,n_j = dt$n,d = d,
								   	 L = L, K = K, J = dt$J,
								     mu_0 = base_measure$mu_0, 
								     kappa_0 = base_measure$kappa_0,
									 nu_0 = base_measure$nu_0, 
									 sigma_0 = base_measure$sigma_0,
									 alpha = a_alpha,rho = a_rho,
									 iter_max = iter_max, warm_up = burn_in,
									 thin = thin, seed = seed, chain = 1,
									 num_posterior_samples = num_posterior_samples)
		}
	}else{
		fit <- beta_nd_nhpp_fit(r= r,
								n_j = dt$n, d = d,
								L = L, K = K, J = dt$J,
								a_0 = base_measure$a_0,
								b_0 = base_measure$b_0,
								mu_sd = base_measure$mu_sd,
								tau_sd = base_measure$tau_sd,
								a_alpha = a_alpha, b_alpha = b_alpha,
								a_rho = a_rho, b_rho = b_rho,
								iter_max = iter_max, warm_up = burn_in,
								thin = thin, seed = seed, chain = 1L,
								num_posterior_samples)

	}


	obj <- list(fit = fit,
				base_measure = base_measure,
				J = dt$J,
				d = if(base_measure$measure=="normal") dt$R*pnorm(d) else dt$R *d,
				n = sum(tail(dt$n,1)),
				call = call,
				R = dt$R,
				r = dt$r,
				K = K,
				L = L)

    out <- ndp(obj)
}

