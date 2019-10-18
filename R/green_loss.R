#' @export
green_loss <- function(object, tau = 0.5)
	UseMethod("green_loss")


#' Lau and Green posterior loss function
#'
#'
#" @param object an hmm object
#' @param tau loss function tuning parameter
#' @result list of two objects (1) the partition corresponding to the minimal loss and
#' (2) the loss value itself
#'
#'@export
green_loss.hmm <- function(object, tau = 0.5){

	loss <- green_loss_engine(object$cluster_assignment[[1]],object$pmat,tau)

	ix <- which.max(loss)

	out <- list(loss = max(loss), mode = object$cluster_assignment[[1]][ix,] )
}
