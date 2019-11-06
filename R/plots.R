#' Plot method for ndp objects
#'
#' The \code{plot} method for ndp-objects provides convenient
#' access to several types of plots useful for visualizing
#' ndp models using the default settings.
#'
#' @export
#' @param  x an ndp object
#' @param plotfun one of "cluster","global", or "network"
#' @param ... ignored
#' @importFrom ggplot2 ggplot aes_string xlab %+replace% theme
#'
plot.ndp <- function(x, plotfun = "cluster",...) {

	if(plotfun == "cluster")
		plot_cluster_densities(x)
	else{
		stop("Not implemented yet")
	}
}

#' plots pairwise probability clustering plot
#' 
#' @template rodriguez
#'
#' @export
#' @param x ndp object
#' @param sort boolean asking whether sorting algorithm should be used to sort
#' pairwise probablity
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn
#' @return ggplot plot object
#' @seealso the supplementary section of the reference for the sorting algorithm.
#'
plot_pairs <- function(x,sort = FALSE)
    UseMethod("plot_pairs")

#' Plots cluster intensity function densities
#'
#' @export
#' @param x ndp object
#' @param p the probability for the credible interval
#' @param pi_threshold intensities with probability of assignment greater than pi_threshold are plotted
#' default is zero
#' @param switch one of two character values "facet" or "color" denoting whether to separate clusters by facet or by color
#' @param transform boolean denoting whether or not to transform distances
#' @return ggplot plot object
#'
plot_cluster_densities <- function(x, p = .9, pi_threshold = .1, switch = "facet",transform = TRUE)
    UseMethod("plot_cluster_densities")

#' Plots global density function 
#'
#' @export
#' @param x ndp object
#' @param p probability mass contained in uncertainty interval
#' @param r optional vector of observed distances
#' @param transform boolean denoting whether or not to transform distances
#' @return ggplot plot object
#' 
plot_global_density <- function(x, p = .9,r = NULL, transform = TRUE)
	UseMethod("plot_global_density")

#' Traceplots of various NDP-NHPP parameters
#'
#' @param x ndp object
#' @param par character vector of parameter names, defaults to c("alpha")
#' @return ggplot plot object
#' @export
plot_traceplots <- function(x, par="alpha")
    UseMethod("plot_traceplots")

#' Network Cluster Plot
#'
#' @export
#' @param x ndp object
#' @param mode_label boolean value to indicate that the network graph should be colored with the Mode cluster label assignments
#' @return ggplot plot object
#' @export
plot_network <- function(x,mode_label = TRUE)
  UseMethod("plot_network")


#' plots pairwise probability clustering plot
#'
#' @export
#' @param x ndp object
#' @param sort boolean asking whether sorting algorithm should be used to sort
#' @return ggplot plot object
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn
#' @importFrom stringr str_replace
#'
plot_pairs.ndp <- function(x,sort = FALSE){

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
#'
#' @export
#' @param x ndp object
#' @param p the probability for the credible interval
#' @param pi_threshold intensities with probability of assignment greater than pi_threshold are plotted
#' default is zero
#' @param switch one of two character values "facet" or "color" denoting whether to separate clusters by facet or by color
#' @param transform boolean value denoting whether or not to transform the distances: FALSE indicates values
#' are shown on the [0,R) domain, while TRUE indicates they are shown on the real line.
#' @return ggplot plot object
#'
plot_cluster_densities.ndp <- function(x, p = .9, pi_threshold = .1, switch = "facet", transform=FALSE){


	### To pass R CMD CHECK
	median <- Intensity_Function <- Distance <- NULL
	Chain <- `Intensity Function` <- quantile <- NULL
	Density <- med <- lower <- upper <- NULL
	###

    ks_to_keep <- which(apply(as.matrix(x$pi),2,median)>pi_threshold)
	if(transform)
		xlabel  <- "Distance (Transformed)"
	else
		xlabel  <- "Distance"


    plt <- x$if_df %>%
		dplyr::mutate(Intensity_Function = factor(Intensity_Function),
					  Distance = (transform==FALSE)*Distance + (transform==TRUE)* stats::qnorm(Distance / ceiling(max(x$if_df$Distance)))) %>%
		dplyr::filter(Intensity_Function %in% ks_to_keep) %>%
		dplyr::rename(`Intensity Function` = Intensity_Function) %>%
		dplyr::group_by(Chain,`Intensity Function`,Distance) %>%
        dplyr::summarise(lower = stats::quantile(Density,.5 - p/2,na.rm=T),
                         med = stats::median(Density,na.rm=T),
                         upper = stats::quantile(Density,.5 + p /2,na.rm=T)) %>%
    ggplot(aes(x=Distance,y=med)) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_blank())

	if(switch == "color"){
		plt <- plt + ggplot2::geom_line(aes(color=`Intensity Function`)) +
			ggplot2::labs(title = "Unnormalized Cluster Intensity Functions", y = "Density", x = xlabel)
	}
	else{
		plt <-  plt + ggplot2::geom_line() + ggplot2::geom_ribbon(aes(ymin= lower,ymax=upper),alpha=0.3) +
			ggplot2::facet_wrap( ~ `Intensity Function`) +
			ggplot2::labs(title = "Normalized Cluster Intensity Functions",
						  subtitle = paste0("Shaded area indicates ",p * 100,"% Credible Interval"),
						  y = "Density", x = xlabel)
	}

    return(plt)
}


#' Traceplots of various NDP-NHPP parameters
#'
#' @export
#' @param x ndp object
#' @param par character vector of parameter names, defaults to c("alpha")
#' @return ggplot plot object
#' @importFrom ggmcmc ggs_traceplot ggs
#'
plot_traceplots.ndp <- function(x,par="alpha"){

    p <- ggs_traceplot( ggs( x[[par]] ) )

    return(p)
}


#' Network Cluster Plot
#'
#' @export
#' @param x ndp object
#' @param mode_label boolean value to indicate that the network graph should be colored with the Mode cluster label assignments
#' @return ggplot plot object
#'
plot_network.ndp <- function(x,mode_label = TRUE){

	### To pass R CMD Check
	nodes <- weight <- `Mode Labels` <- NULL
	###


    t <- tidygraph::as_tbl_graph(x$pmat,directed=FALSE)
	if(mode_label){
		mode <- assign_mode(x)
		if(any(mode==0))
			mode  <- mode + 1
		t <- t %>% tidygraph::activate(nodes) %>% tidygraph::mutate(`Mode Labels` = factor(mode) )
	}
    p <- t %>% ggraph::ggraph(layout='nicely') +
        ggraph::geom_edge_link(aes(alpha=weight)) +
		 ggplot2::ggtitle("Network Plot") + ggplot2::theme_void()
	if(mode_label)
		p <- p + ggraph::geom_node_point(aes(color=`Mode Labels`))
	else
		p <- p + ggraph::geom_node_point()

    return(p)
}

#' Plots global density function 
#'
#' Plots global density function 
#'
#' @export
#' @param x ndp object
#' @param p probability mass contained in uncertainty interval
#' @param r optional vector of observed distances
#' @param transform boolean denoting whether or not to transform distances
#' @return ggplot plot object
#'
plot_global_density.ndp <- function(x, p = 0.9, r = NULL,transform = FALSE ){

	### to pass R CMD CHECK
	med <- lower <- upper <- Global_Density <- NULL
	distances <- Chain <- Distance <- NULL
	###

	if(p >= 1 || p <= 0 )
		stop("p must be in (0,1)")

    plt <- x$global_density %>%
		dplyr::mutate(Distance = (transform==FALSE)*Distance + (transform==TRUE)*stats::qnorm(Distance / ceiling(max(x$if_df$Distance)))) %>%
		dplyr::group_by(Chain,Distance) %>%
        dplyr::summarise(lower = stats::quantile(Global_Density,.5 - p / 2,na.rm=T),
                         med = stats::median(Global_Density,na.rm=T),
                         upper = stats::quantile(Global_Density,.5 + p /2,na.rm=T)) %>%
    ggplot(aes(x=Distance,y=med)) + ggplot2::geom_line() +
    ggplot2::geom_ribbon(aes(ymin= lower,ymax=upper),alpha=0.3) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~Chain) +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    ggplot2::labs(title = "Global Density NDP Estimate",
         subtitle = paste0("Shaded area indicates ",p*100,"% Credible Interval"),
         y = "Density")
	if(is.null(r))
		return(plt)
	else if(is.numeric(r)){
		p2 <- dplyr::tibble(distances = r) %>%  ggplot2::ggplot(aes(x=distances)) +
			ggplot2::geom_density() + ggplot2::theme_bw() +
			ggplot2::labs(title = "Global Kernel Density Estimate")
		return(gridExtra::grid.arrange(plt,p2,nrow=1))
	}
}
