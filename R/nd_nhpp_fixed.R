#j#' Estimate the nonhomgogenous poisson process intensity function from grouped data using a normal kernel
#j#'
#j#' @param distances_col name of column in data that contains distances
#j#' @param id_col name of column in data  that contains id grouping variable
#j#' @param data data.frame object that contains grouped distances
#j#' @param L component truncation number
#j#' @param K intensity cluster truncation number
#j#' @param mu_0 mean hyperparameter for mu normal base measure. Default is 0; Normal(0,1).
#j#' @param kappa_0 variance hyperparameter for mu normal base measure. Default is 1; Normal(0,1).
#j#' @param nu_0 df hyperparameter for sigma inv chisq base measure. Default is 1; InvChisq(1,1);
#j#' @param sigma_0 scale hyperparameter for sigma inv chisq base measure. Default is 1; InvChisq(1,1);
#j#' @param alpha outer DP concentration parameter
#j#' @param rho inner DP concentration parameter
#j#' @param iter_max total number of iterations for which to run sampler
#j#' @param warm_up number of iterations for which to burn-in or "warm-up" sampler
#j#' @param thin number of iterations to thin by
#j#' @param seed integer with which to initialize random number generator
#j#'
#j#' @export
#jnd_nhpp_fixed <- function(distances_col,id_col,
#j						  data = NULL,
#j						  mu_0 = 0, kappa_0 = 1,
#j						  nu_0 = 1, sigma_0 = 1,
#j						  alpha = 1, rho = 1,
#j						  L = 4 , K = 4,
#j						  iter_max, warm_up,
#j						  thin = 1,
#j						  seed = NULL) {
#j
#j    call <- match.call(expand.dots=TRUE)
#j	r <- data %>% dplyr::arrange(!! dplyr::sym(id_col)) %>%
#j		dplyr::select(!!!distances_col) %>% dplyr::pull()
#j	n_j <- data %>%  dplyr::arrange(!! dplyr::sym(id_col)) %>%
#j		dplyr::group_by(!! dplyr::sym(id_col)) %>% dplyr::count() %>%
#j		dplyr::ungroup() %>% dplyr::mutate(start = cumsum(n)) %>%
#j		dplyr::mutate(start_ = tidyr::replace_na(dplyr::lag(start),0)) %>%
#j		dplyr::select(-start) %>% dplyr::rename(start = start_,go = n) %>%
#j		dplyr::select(start,go) %>% as.matrix()
#j    J <-  nrow(n_j)
#j    ## basic checks
#j	stopifnot((iter_max > warm_up) && (warm_up > 0))
#j	stopifnot(L>0 && K >0)
#j    if(any(r<=0))
#j        stop("all r must be positive numbers", .call = FALSE)
#j    if(any(r>=1)){
#j        R <- ceiling(max(r)) ## maximum radius
#j        r_ <- r / R ## scaled event radii
#j    }else{
#j        r_ <-  r
#j        R <- 1
#j    }
#j
#j
#j
#j    if(is.null(seed))
#j        seed <- 1L
#j
#j	r_  <-  stats::qnorm(r_)
#j
#j    d <- seq(from = floor(min(r_)), to = ceiling(max(r_)), by = 0.01) ## distance grid
#j    num_posterior_samples <- sum(seq(from=warm_up+1,to = iter_max,by=1) %% thin == 0 )
#j    fit <- list(nd_nhpp_fixed_fit(r= r_, n_j = n_j, d = d,
#j								  L = L, K = K, J = J,
#j								  mu_0 = mu_0, kappa_0 = kappa_0,
#j								  nu_0 = nu_0, sigma_0 = sigma_0,
#j								  alpha = alpha, rho = rho,
#j								  iter_max = iter_max, warm_up = warm_up,
#j								  thin = thin, seed = seed, chain = 1,
#j								  num_posterior_samples = num_posterior_samples))
#j    d <- stats::pnorm(d)
#j
#j    out <- ndp_fixed(c(list(K = K, L = L, d = R*d,
#j                      n = sum(n_j[,2]), call = call),fit),1,alpha,rho,J)
#j}
#j
