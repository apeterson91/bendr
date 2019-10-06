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
assign_mode.ndp <- function(x,ics = NULL){


	error <- get_error(x,ics)
	if(!is.null(ics))
		return(x$cluster_assignment[[1]][which.min(error),ics])
	else
		return(x$cluster_assignment[[1]][which.min(error),])
}

#' Calculates Error Distribution
#'
#' @export
#' @method get_error ndp
#' @param x ndp object
#' @return vector of error corresponding to specific iteration's cluster configuration
#'
get_error <- function(x,ics = NULL)
	UseMethod("get_error")


#' @export
#'
get_error.ndp <- function(x, ics = NULL){

	A <- x$pmat

	A[upper.tri(A)] <- A[lower.tri(A)]

	if(!is.null(ics)){
		A <- A[ics,ics]
		error <- sapply(1:nrow(x$cluster_assignment[[1]]),function(y) sum((A - get_adj_mat(x$cluster_assignment[[1]][y,ics]))^2)  )
	}
	else
		error <- sapply(1:nrow(x$cluster_assignment[[1]]),function(y) sum((A - get_adj_mat(x$cluster_assignment[[1]][y,]))^2)  )

	return(error)
}


#--- internal only

get_adj_mat <- function(x){

	A <- matrix(nrow = length(x),ncol = length(x))
	for(i in 1:length(x)){
		for(j in 1:length(x))
			A[,i] <- (x == x[j])*1
	}
	return(A)
}


