#' Plot method for ndp objects
#'
#' The \code{plot} method for ndp-objects provides convenient
#' access to several types of plots useful for visualizing
#' ndp models using the default settings.
#'
#' @export
#' @param  x an ndp object
#' @param plotfun one of c("cluster","pairs","global","trace" )
#' @param ... optional arguments for plot function 
#' @seealso \code{\link{plot_cluster_densities}} \code{\link{plot_pairs}} \code{\link{plot_global_density}} \code{\link{traceplot}}
#'
plot.ndp <- function(x, plotfun = "cluster",...) {

	stopifnot(plotfun %in% c("cluster","pairs","global","trace" ))

	p <- switch(plotfun,
		   "cluster" = plot_cluster_densities(x,...),
		   "pairs" = plot_pairs(x,...),
		   "global" = plot_global_density(x,...),
		   "trace" = traceplot(x,...)
	   )
	return(p)

}
