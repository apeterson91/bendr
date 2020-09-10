#' ndp objects
#'
#' @param object A list provided by the \code{\link{bend}} function
#' @return An ndp object for easy model interrogation
#'
ndp <- function(object){


	K <- object$K
	pmat <- object$fit$cluster_pair_probability
	colnames(pmat) <- paste0("Group_",1:ncol(pmat))
	p <- object$fit$pi_samples
	w <- object$fit$w_samples
	mu <- object$fit$mu_samples
	tau <- object$fit$tau_samples
	alpha <- object$fit$alpha_samples
	rho <- object$fit$rho_samples
	dgd <- object$fit$global_density
	dkd <- object$fit$intensities
	colnames(dgd) <- as.character(1:length(object$d))
	gdensity <- purrr::map_dfr(1:nrow(dgd),function(s){
								   dplyr::tibble(sim_id = s,
												 Distance = object$d,
												  Density = dgd[s,])
							})
	dlen <- length(object$d)
	kdensity <- purrr::map_dfr(0:(K-1),function(k) { 
					purrr::map_dfr(1:nrow(dgd),function(s){ 
								   dplyr::tibble(sim_id = s,
												 Distance = object$d,
												 K = k+1,
												 Density = dkd[s,(k*dlen+1):(k*dlen+dlen) ] ) })
				})

	gd <- expand.grid(K = 1:K,
					  L=1:object$L)
	gdnms <- paste0("K: ", gd$K,"L: ",gd$L)
	paste_nms <- function(nm){
		paste0(nm, " ",gdnms)
	}
	colnames(p) <- paste0("pi K:", 1:K)
	colnames(w) <- paste_nms("w")
	colnames(mu) <- paste_nms("mu")
	if(is_normal(object$base_measure))
		colnames(tau) <- paste_nms("tau")
	else
		colnames(tau) <- "tau"
	colnames(alpha) <- "alpha"
	colnames(rho) <- "rho"

	out <- list(cluster_assignment = object$fit$cluster_assignment,
				component_assignment = object$fit$component_assignment,
				pmat = pmat,
				call = object$call,
				base_measure = object$base_measure,
				w = w,
				p = p,
				mu = mu,
				tau = tau,
				alpha = alpha,
				rho = rho,
				gdensity = gdensity,
				kdensity = kdensity,
				n = object$n,
				K = object$K,
				J = object$J,
				r = object$r,
				R = object$R
				)

    structure(out, class = c("ndp"))
}
