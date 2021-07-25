#' Methods for ndp objects
#' 
#' Common methods that allow for easy model interrogation
#' 
#' @name ndp-methods
#' @aliases nsamples 
#' 
#' @param object ndp object
#' @param ... Ignored
#' 

#'
#' @rdname ndp-methods
#' @export
#' @importFrom rstantools nsamples
nsamples.ndp <- function(object, ...) {
  nrow(object$alpha)
}


#' Cluster Credible Ball
#'
#' Calculates the credible ball (analog to the credible interval)  for the cluster assignments
#' @export 
#' @param object an ndp object
#' @seealso Wade, Sara, and Zoubin Ghahramani. "Bayesian cluster analysis: Point estimation and credible balls (with discussion)." Bayesian Analysis 13.2 (2018): 559-626.
#' 
clustercb <- function(object) UseMethod("clustercb")


#' @describeIn consensus
#' @export
clustercb.ndp <- function(object){


	cluster_assignment <- object$cluster_assignment +  1

	P <- makeSymm(object$pmat)

	diag(P) <- 1

	NDP.VI <- mcclust.ext::minVI(P,cluster_assignment,
							method=("all"),
							include.greedy = TRUE,
							suppress.comment = FALSE)
	ndp.cb <- mcclust.ext::credibleball(NDP.VI$cl[4,],cluster_assignment)

	return(list(CredibleBall = ndp.cb, VI = NDP.VI))

}




# Internal ------------------------------------

makeSymm <- function(m) {

	  m[upper.tri(m)] <- t(m)[upper.tri(m)]
	  return(m)

}
