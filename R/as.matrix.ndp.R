#' Retrieve Parameter Samples in Matrix Form
#'
#' Returns the parameter samples from the model associated with the 
#' means, variances, concentration parameters and NDP probabilities
#' as  a S X P matrix where S is the number of samples and P is the total number of parameters
#' 
#' @param x a ndp object
#' @param ... ignored
#' @export
#' @seealso as.array.stapDP for an array like container for the NDP components
#'
as.matrix.ndp <- function(x,...){

	return(cbind(x$p,x$w,x$mu,x$tau,x$alpha,x$rho))
}

#' Retrieve NDP parameter samples in Array Form
#'
#' @param x a ndp object
#' @param ... ignored
#' @export
#' @seealso as.matrix.ndp for a matrix container of *all* parameters
#' 
as.array.ndp <- function(x,...){

	mu <- abind::abind(lapply(1:x$K,function(k) x$mu[,grep(x = colnames(mu),paste("K:",k),value=T)]),along=3)
}
