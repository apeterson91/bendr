#' Lau and Green posterior loss function
#'
#' @param object an  object with N X N cluster pairwise probability matrix and N X K cluster matrix list entries
#' @param truth True adjacency matrix
#' @param tau loss function tuning parameter
#' @param a misclassification penalty parameter
#' @param b misclassification penalty parameter - see Green and Lau for details
#' @return list of three objects (1)loss:  the loss corresponding to each clustering 
#' configuration, (2) best_loss_ix:  The index corresponding to the "best" partition,
#' and (3) mode: the best partition itself
#' @export
green_loss <- function(object, truth = NULL, tau = 0.5, a = 1, b = 1)
	UseMethod("green_loss")


#' @describeIn green_loss classification loss function
#' @export
green_loss.default <- function(object, truth = NULL, tau = 0.5, a = 1, b = 1){

	if(is.null(truth)){
		loss <- green_loss_unknown(object$cluster_assignment,object$pmat,tau)
		ix <- which.max(loss)
	}
	else{
		stopifnot(dim(truth)[1] == ncol(object$cluster_assignment))
		stopifnot(dim(truth)[2] == ncol(object$cluster_assignment))
		loss <- green_loss_known(object$cluster_assignment,object$pmat,truth,a,b)
		ix <- which.min(loss)
	}

	out <- list(loss = loss,best_loss_ix = ix , mode = object$cluster_assignment[ix,] )
}

