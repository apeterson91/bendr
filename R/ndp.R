#' ndp objects
#'
#' @param object A list provided by the \code{\link{bend}} function
#' @return An ndp object for easy model interrogation
#'
ndp <- function(object){


	K <- object$K
	p <- object$fit$pi_samples
	w <- object$fit$w_samples
	mu <- object$fit$mu_samples
	tau <- object$fit$tau_samples
	alpha <- object$fit$alpha_samples
	rho <- object$fit$rho_samples
	tau <- object$fit$tau
	gd <- expand.grid(K = 1:K,L=1:object$L)
	gdnms <- paste0("K: ", gd$K,"L: ",gd$L)
	paste_nms <- function(nm){
		paste0(nm,gdnms)
	}
	colnames(p) <- paste0("pi K:", 1:K)
	colnames(w) <- paste_nms("w")
	colnames(mu) <- paste_nms("mu")
	colnames(tau) <- paste_nms("tau")
	colnames(alpha) <- "alpha"
	colnames(rho) <- "rho"



	out <- list(cluster_assignment = object$fit$cluster_assignment,
				component_assignment = object$fit$component_assignment,
				pmat = object$fit$cluster_pair_probability,
				w = w,
				p = p,
				mu = object$fit$mu_samples,
				tau = object$fit$tau_samples,
				alpha = object$fit$alpha_samples,
				rho = object$fit$rho_samples,
				alpha_prior = object$fit$alpha_prior,
				rho_prior = object$fit$rho_prior,
				gdensity = object$fit$global_intensity,
				cdensity = object$fit$intensities
				)

    structure(out, class = c("ndp"))
}
