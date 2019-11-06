#' Plot method for bndp objects
#'
#' The \code{plot} method for bndp-objects provides convenient
#' access to several types of plots useful for visualizing
#' bndp models using the default settings.
#'
#' @export
#' @param  x an bndp object
#' @param plotfun one of "cluster","global", or "network"
#' @param ... ignored
#' @importFrom ggplot2 ggplot aes_string xlab %+replace% theme
#'
plot.bndp <- function(x, plotfun = "cluster",...) {

	if(plotfun == "cluster")
		plot_cluster_densities(x)
	else{
		stop("Not implemented yet")
	}
}


#' plots pairwise probability clustering plot
#' 
#' @export
#' @param x bndp object
#' @param sort boolean asking whether sorting algorithm should be used to sort
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn
#' @return ggplot plot object
#' 
plot_pairs.bndp <- function(x,sort = FALSE){

	### To pass R CMD Check
	Group_2 <- Group_1 <- Probability <- NULL
	###

	P <- x$pmat

	makeSymm <- function(m) {
	  m[upper.tri(m)] <- t(m)[upper.tri(m)]
	  return(m)
	}

	P <- makeSymm(P)
	if(sort){

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
	}
	else{
		A <- 1:nrow(P)
	}


	p <- dplyr::as_tibble(P[A,A]) %>% dplyr::mutate(Group_1 = 1:dplyr::n()) %>%
	  tidyr::gather(dplyr::contains("V"), key = "Group_2", value = "Probability") %>%
	  dplyr::mutate("Group_2" = as.numeric(stringr::str_replace(Group_2,"V",""))) %>%
	  ggplot(aes(x=Group_1,y=Group_2,fill=Probability)) +
	  geom_tile() + scale_fill_gradientn(colours=c("white","grey","black"),limits=c(0,1)) +
	  ggplot2::theme_bw() + ggplot2::labs(title = "Pairwise Probability of Function Clustering",
										  x="Group 1", y = "Group 2") +
	  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
			panel.grid.minor = ggplot2::element_blank())

    return(p)
}

#' Plots cluster intensity function densities
#' @export
#' @param x bndp object
#' @param p the probability for the credible interval
#' @param pi_threshold intensities with probability of assignment greater than pi_threshold are plotted
#' default is zero
#' @param switch one of two character values "facet" or "color" denoting whether to separate clusters by facet or by color
#' @param transform not used for objects of type bndp
#' @return ggplot plot object
#'
plot_cluster_densities.bndp <- function(x, p = .9,pi_threshold = .1, switch = "facet", transform = FALSE){

	### To pass R CMD CHECK
	median <- Intensity_Function <- Distance <- NULL
	Chain <- `Intensity Function` <- quantile <- NULL
	Density <- med <- lower <- upper <- transform <- NULL
	###

    ks_to_keep <- which(apply(as.matrix(x$pi),2,median)>pi_threshold)
	xlabel  <- "Distance"


    p <- x$if_df %>% 
		dplyr::mutate(Intensity_Function = factor(Intensity_Function),
					  Distance = Distance)   %>% 
		dplyr::filter(Intensity_Function %in% ks_to_keep) %>%
		dplyr::rename(`Intensity Function` = Intensity_Function) %>% 
		dplyr::group_by(Chain,`Intensity Function`,Distance) %>%
        dplyr::summarise(lower = quantile(Density,.5 + p/2,na.rm=T),
                         med = median(Density,na.rm=T),
                         upper = quantile(Density,.5 + p /2,na.rm=T)) %>%
    ggplot(aes(x=Distance,y=med)) + 
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_blank())
    
	if(switch == "color"){
		p <- p + ggplot2::geom_line(aes(color=`Intensity Function`)) + 
			ggplot2::labs(title = "Unormalized Cluster Intensity Functions", y = "Density", x = xlabel)
	}
	else{
		p <-  p + ggplot2::geom_line() + ggplot2::geom_ribbon(aes(ymin= lower,ymax=upper),alpha=0.3) +
			ggplot2::facet_wrap( ~ `Intensity Function`) + 
			ggplot2::labs(title = "Unnormalized Cluster Intensity Functions",
						  subtitle = paste0("Shaded area indicates ",p*100,"% Credible Interval"),
						  y = "Density", x = xlabel)
	}
    
    return(p)
}


#' Network Cluster Plot
#'
#' @export 
#' @param x bndp object
#' @param mode_label boolean value to indicate that the network graph should be colored with the Mode cluster label assignments
#' @return ggplot plot object
#'
plot_network.bndp <- function(x, mode_label = TRUE){
	plot_network.ndp(x, mode_label = mode_label)
}


#' Plots global density function 
#'
#' @export
#' @param x bndp object
#' @param p probability mass contained in uncertainty interval
#' @param r optional vector of observed distances
#' @param transform boolean denoting whether or not to transform distances
#' @return ggplot plot object
#'
plot_global_density.bndp <- function(x, p = 0.9, r = NULL, transform = TRUE){


	### to pass R CMD CHECK
	med <- lower <- upper <- Global_Density <- NULL
	distances <- Chain <- Distance <- NULL
	###

	if(p >= 1 || p <= 0 )
		stop("p must be in (0,1)")

    p <- x$global_density %>% dplyr::group_by(Chain,Distance) %>%
        dplyr::summarise(lower = quantile(Global_Density,.5 - p / 2,na.rm=T),
                         med = median(Global_Density,na.rm=T),
                         upper = quantile(Global_Density,.5 + p /2,na.rm=T)) %>%
    ggplot(aes(x=Distance,y=med)) + ggplot2::geom_line() +
    ggplot2::geom_ribbon(aes(ymin= lower,ymax=upper),alpha=0.3) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~Chain) +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    ggplot2::labs(title = "Global Density NDP Estimate",
         subtitle = paste0("Shaded area indicates ",p,"% Credible Interval"),
         y = "Density")
	if(is.null(r))
		return(p)
	else if(is.numeric(r)){
		p2 <- dplyr::tibble(distances = r) %>%  ggplot2::ggplot(aes(x=distances)) +
			ggplot2::geom_density() + ggplot2::theme_bw() +
			ggplot2::labs(title = "Global Kernel Density Estimate")
		return(gridExtra::grid.arrange(p,p2,nrow=1))
	}
}
