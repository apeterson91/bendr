#' Estimates Posterior Mode Cluster Assignment
#'
#' @export
#' @method assign_mode ndp
#' @param x ndp object
#' @return vector of cluster assignments
#'
assign_mode <- function(x,ics = NULL)
	UseMethod("assign_mode")


#' @export
#'
assign_mode.ndp <- function(x){


	error <- get_square_error(x)
	return(x$cluster_assignment[[1]][which.max(error)])
}

#' Calculates Error Distribution under square loss function
#'
#' @export
#' @method get_square_error ndp
#' @param x ndp object
#' @return vector of error corresponding to specific iteration's cluster configuration
#'
get_square_error <- function(x)
	UseMethod("get_square_error")


#' @export
#'
get_square_error.ndp <- function(x){

	error <- square_error(x$cluster_assignment[[1]],x$pmat)

	return(error)
}


