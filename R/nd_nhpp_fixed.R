#' Estimate the nonhomgogenous poisson process intensity function from grouped data using a normal kernel
#'
#' @param distances_col name of column in data that contains distances
#' @param id_col name of column in data  that contains id grouping variable
#' @param data data.frame object that contains grouped distances
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
nd_nhpp_fixed <- function(distances_col,id_col,
						  data = NULL,
						  mu_0 = 0, kappa_0 = 1,
						  nu_0 = 1, sigma_0 = 1,
						  alpha = 1, rho = 1,
						  L = 4 , K = 4,
						  iter_max, warm_up,
						  thin = 1,
						  seed = NULL) {

    call <- match.call(expand.dots=TRUE)
	r <- data %>% dplyr::arrange(!! dplyr::sym(id_col)) %>%
		select(!!!distances_col) %>% pull()
	n_j <- data %>%  dplyr::arrange(!! dplyr::sym(id_col)) %>%
		dplyr::group_by(!! dplyr::sym(id_col)) %>% dplyr::count() %>%
		dplyr::ungroup() %>% dplyr::mutate(start = cumsum(n)) %>%
		mutate(start_ = replace_na(dplyr::lag(start),0)) %>%
		dplyr::select(-start) %>% dplyr::rename(start = start_,go = n) %>%
		dplyr::select(start,go) %>% as.matrix()
    J <-  nrow(n_j)
    ## basic checks
	stopifnot((iter_max > warm_up) && (warm_up > 0))
	stopifnot(is.integer(iter_max) && is.integer(warm_up))
	stopifnot(L>0 && K >0)
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

