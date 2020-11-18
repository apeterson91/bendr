
#' plots pairwise probability clustering plot
#'
#' @template rodriguez
#'
#' @export
#' @param x object with upper triangle pairwise probability matrix list entry
#' @param sample positive integer denoting size of random sample for which the pairwise graph can be subset - useful for speeding up graphing
#' @param sort boolean asking whether sorting algorithm should be used to sort
#' pairwise probablity
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn theme xlab ylab
#' @importFrom stringr str_replace
#' @return ggplot plot object
#' @seealso the supplementary section of the reference for the sorting algorithm.
#'
plot_pairs <- function(x, sample = NULL, sort = FALSE)
    UseMethod("plot_pairs")

#' Plots cluster intensity function densities
#'
#' @export
#' @param x ndp object
#' @param p the probability for the credible interval
#' @param pi_threshold intensities with probability of assignment greater than pi_threshold are plotted
#' default is zero
#' @param style one of two character values "facet" or "color" denoting whether to separate clusters by facet or by color
#' @param transform boolean denoting whether or not to transform distances
#' @return ggplot plot object
#'
plot_cluster_densities <- function(x, p = .9, pi_threshold = .1, style = "facet",transform = TRUE)
    UseMethod("plot_cluster_densities")

#' Plots global density function
#'
#' @export
#' @param x ndp object
#' @param p probability mass contained in uncertainty interval
#' @param transform boolean denoting whether or not to transform distances
#' @return ggplot plot object
#'
plot_global_density <- function(x, p = .9, transform = TRUE )
	UseMethod("plot_global_density")

#' Traceplots of various NDP-NHPP parameters
#'
#' @param x ndp object
#' @param par character vector of parameter names, defaults to c("alpha")
#' @return ggplot plot object
#' @export
traceplot <- function(x, par="alpha")
    UseMethod("traceplot")


#' plots pairwise probability clustering plot
#'
#' @describeIn plot_pairs constructs pairwise probability heat map
#' @export
#'
plot_pairs.default <- function(x,sample = NULL, sort = FALSE){

	### To pass R CMD Check
	Group_2 <- Group_1 <- index <- Probability <- NULL
	###
	if(is.null(sample))
		P <- x$pmat
	else if(is.integer(sample) && sample>0){
		ics <- sample(1:nrow(x$pmat),sample)
		P <- x$pmat[ics,ics]
	}else
		stop("sample is not an integer")

	makeSymm <- function(m) {
	  m[upper.tri(m)] <- t(m)[upper.tri(m)]
	  return(m)
	}

	P <- makeSymm(P)
	if(sort){
	    colnames(P) <- NULL
	    rownames(P) <- NULL

		i <- which(P == max(P),arr.ind = T)[1]
		j <- which(P == max(P),arr.ind = T)[2]

		A <- c(i,j)
		Omega <- 1:nrow(P)
		A_c <- setdiff(Omega,A)
		while(length(A_c)){
		  probs <- sapply(A_c,function(x) max(P[x,A]))
		  j <- A_c[which.max(probs)]
		  A <- c(A,j)
		  A_c <- setdiff(Omega,A)
		}
		P <- P[A,A]
	}
	else{
		A <- 1:nrow(P)
		P <- P[A,A]
	}
	colnames(P) <- paste0("ID_",1:ncol(P))


	p <- suppressMessages(dplyr::as_tibble(P)) %>%
		dplyr::mutate(Group_1 = 1:dplyr::n()) %>%
		tidyr::pivot_longer(dplyr::contains("ID_"), names_to = "Group_2", values_to = "Probability") %>%
		dplyr::mutate("Group_2" = as.numeric(stringr::str_replace(Group_2,"ID_",""))) %>%
		ggplot(aes(x=Group_1,y=Group_2,fill=Probability)) +
		geom_tile() + scale_fill_gradientn(colours=c("white","grey","black"),limits=c(0,1)) +
		ggplot2::theme_bw() + ggplot2::labs(title = "Pairwise Probability of Function Clustering",
										  x="Subject 1", y = "Subject 2") +
		ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
					   panel.grid.minor = ggplot2::element_blank())


    return(p)
}

#'
#' @export
#' @describeIn plot_cluster_densities cluster density estimates
#'
plot_cluster_densities.ndp <- function(x, p = .9, pi_threshold = .1, style = "facet", transform = FALSE){


	### To pass R CMD CHECK
	 K <- Distance <- NULL
	Density <- Median <- Lower <- Upper <- NULL
	###

	stopifnot(p>=0 && p<=1)
	l <-  .5 - p/2
	u <- .5 + p/2

    ks_to_keep <- which(apply(x$p,2,median)>pi_threshold)
	tfun <- function(y) y
	xlabel  <- "Distance"
	if(transform){
		if(x$base_measure$measure == "Normal")
			xlabel  <- "Distance (Transformed)"
			tfun <- function(y){
				qnorm(y/x$R)
			}
	}

	x$kdensity %>%
		dplyr::filter(K %in% ks_to_keep) %>%
		dplyr::mutate(Distance = tfun(Distance),
					  K = factor(K)) %>%
		dplyr::group_by(K,Distance) %>%
		dplyr::summarise(Lower = quantile(Density,l,na.rm=TRUE),
						 Median = quantile(Density,.5,na.rm=TRUE),
						 Upper = quantile(Density,u,na.rm=TRUE))-> pltdf

    pltdf %>% ggplot(aes(x=Distance,y=Median)) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_blank()) -> plt

	if(style == "color"){
		plt <- plt + ggplot2::geom_line(aes(color=K )) +
			ggplot2::labs(title = "Unnormalized Cluster Intensity Functions",
						  y = "Density",
						  x = xlabel)
	}
	else{
		plt <-  plt + ggplot2::geom_line() +
			ggplot2::geom_ribbon(aes(ymin= Lower,ymax=Upper),alpha=0.3) +
			ggplot2::facet_wrap( ~ K) +
			ggplot2::labs(title = "Normalized Cluster Intensity Functions",
						  subtitle = paste0("Shaded area indicates ",p * 100,"% Credible Interval"),
						  y = "Density", x = xlabel)
	}

    return(plt)
}


#'
#' @export
#' @describeIn traceplot parameter traceplots
#'
traceplot.ndp <- function(x,par="alpha"){

	mat <- as.matrix(x)
	if(!is.null(par)){
		stopifnot(par %in% colnames(mat))
	}
	mat <- mat [,par,drop=F]
    p <- bayesplot::mcmc_trace(list(mat))
    return(p)
}


#'
#' Plots global density function
#'
#' @export
#' @describeIn plot_global_density plots global density estimate
#'
#'
plot_global_density.ndp <- function(x, p = 0.9, transform = FALSE){

	### to pass R CMD CHECK
	med <- lower <- upper <- Density <- NULL
	distances <- Chain <- Distance <- NULL
	###

	if(p >= 1 || p <= 0 )
		stop("p must be in (0,1)")
	l <-  .5 - p / 2
	u <- .5 + p / 2

	xlabel  <- "Distance"
	tfun <- function(y) y
	if(transform){
		if(x$base_measure$measure == "Normal")
			xlabel  <- "Distance (Transformed)"
			tfun <- function(y){
				qnorm(y/x$R)
			}
	}

    plt <- x$gdensity %>%
		dplyr::mutate(Distance = tfun(Distance)) %>%
		dplyr::group_by(Distance) %>%
        dplyr::summarise(lower = stats::quantile(Density,l,na.rm=T),
                         med = stats::median(Density,na.rm=T),
                         upper = stats::quantile(Density,u,na.rm=T)) %>%
    ggplot(aes(x=Distance,y=med)) + ggplot2::geom_line() +
    ggplot2::geom_ribbon(aes(ymin= lower,ymax=upper),alpha=0.3) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    ggplot2::labs(title = "Global Density NDP Estimate",
         subtitle = paste0("Shaded area indicates ",p*100,"% Credible Interval"),
         y = "Density")
}
