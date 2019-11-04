#' @export
green_loss <- function(object, truth = NULL, tau = 0.5, a = 1, b = 1)
	UseMethod("green_loss")


#' Lau and Green posterior loss function
#'
#'
#" @param object an hmm object
#' @param truth True adjacency matrix
#' @param tau loss function tuning parameter
#' @param a misclassification penalty parameter
#' @param b misclassification penalty parameter - see Green and Lau for details
#' @return list of two objects (1) the partition corresponding to the minimal loss and
#' (2) the loss value itself
#'
#'@export
green_loss.ndp <- function(object, truth = NULL, tau = 0.5, a = 1, b = 1){

	if(is.null(truth)){
		loss <- green_loss_unknown(object$cluster_assignment[[1]],object$pmat,tau)
		ix <- which.max(loss)
	}
	else{
		stopifnot(dim(truth)[1] == ncol(object$cluster_assignment[[1]]))
		stopifnot(dim(truth)[2] == ncol(object$cluster_assignment[[1]]))
		loss <- green_loss_known(object$cluster_assignment[[1]],object$pmat,truth,a,b)
		ix <- which.min(loss)
	}

	out <- list(loss = loss,best_loss_ix = ix , mode = object$cluster_assignment[[1]][ix,] )
}
