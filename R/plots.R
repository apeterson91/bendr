#' Plot method for ndp objects
#'
#' The \code{plot} method for \link{ndp-objects} provides convenient
#' access to several types of plots useful for visualizing
#' ndp models using the default settings.
#'
#' @method plot ndp
#' @export
#' @param  x an ndp object
#' @param plotfun one of "cluster","global", or "network"
#' @importFrom ggplot2 ggplot aes_string xlab %+replace% theme
#'
plot.ndp <- function(x, plotfun = "cluster") {

	if(plotfun == "cluster")
		plot_cluster_densities(x)
}

#' @export
plot_pairs <- function(x,sort = FALSE)
    UseMethod("plot_pairs")

#' @export
plot_cluster_densities <- function(x, p = .9, pi_threshold = .1, switch = "facet",transform = TRUE)
    UseMethod("plot_cluster_densities")

#' @export
plot_global_density <- function(x, p = .9,r = NULL, transform = TRUE)
	UseMethod("plot_global_density")

#' @export
plot_traceplots <- function(x, par="alpha")
    UseMethod("plot_traceplots")

#' @export
plot_network <- function(x,sample=NULL)
  UseMethod("plot_network")


#' plots pairwise probability clustering plot
#'
#' @export
#' @method plot_pairs ndp
#' @param x ndp object
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn
#' @return ggplot plot object
#'
plot_pairs.ndp <- function(x,sort = FALSE){

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
	  ggplot2::theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank())


    return(p)
}

#' Plots cluster intensity function densities
#' @export
#' @method plot_cluster_densities ndp
#' @param x ndp object
#' @param p the probability for the credible interval
#' @param pi_threshold intensities with probability of assignment greater than pi_threshold are plotted
#' default is zero
#' @param switch one of two character values "facet" or "color" denoting whether to separate clusters by facet or by color
#' @return ggplot plot object
#'
plot_cluster_densities.ndp <- function(x, p = .9,pi_threshold = .1, switch = "facet",transform=TRUE){

    ks_to_keep <- which(apply(as.matrix(x$pi),2,median)>pi_threshold)
	if(transform)
		xlabel  <- "Distance"
	else
		xlabel  <- "Distance (Transformed)"


    plt <- x$if_df %>%
		dplyr::mutate(Intensity_Function = factor(Intensity_Function),
					  Distance = (transform==TRUE)*Distance + (transform==FALSE)*qnorm(Distance / ceiling(max(x$if_df$Distance)))) %>%
		dplyr::filter(Intensity_Function %in% ks_to_keep) %>%
		dplyr::rename(`Intensity Function` = Intensity_Function) %>%
		dplyr::group_by(Chain,`Intensity Function`,Distance) %>%
        dplyr::summarise(lower = quantile(Density,.5 - p/2,na.rm=T),
                         med = median(Density,na.rm=T),
                         upper = quantile(Density,.5 + p /2,na.rm=T)) %>%
    ggplot(aes(x=Distance,y=med)) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_blank())

	if(switch == "color"){
		plt <- plt + ggplot2::geom_line(aes(color=`Intensity Function`)) +
			ggplot2::labs(title = "Cluster Normalized Intensity Functions", y = "Density", x = xlabel)
	}
	else{
		plt <-  plt + ggplot2::geom_line() + ggplot2::geom_ribbon(aes(ymin= lower,ymax=upper),alpha=0.3) +
			ggplot2::facet_wrap( ~ `Intensity Function`) +
			ggplot2::labs(title = "Cluster Normalized Intensity Functions",
						  subtitle = paste0("Shaded area indicates ",p * 100,"% Credible Interval"),
						  y = "Density", x = xlabel)
	}

    return(plt)
}


#' Traceplots of various NDP-NHPP parameters
#'
#' @export
#' @method plot_traceplots ndp
#' @param x ndp object
#' @param par character vector of parameter names, defaults to c("alpha")
#' @return ggplot plot object
#'
plot_traceplots.ndp <- function(x,par="alpha"){

    p <- ggmcmc::ggs_traceplot( ggmcmc::ggs( x[[par]] ) )

    return(p)
}


#' Network Cluster Plot
#'
#' @export
#' @method plot_network ndp
#' @param x ndp object
#' @param mode_label boolean value to indicate that the network graph should be colored with the Mode cluster label assignments
#' @return ggplot plot object
#'
plot_network.ndp <- function(x,mode_label = TRUE){



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

#' Plots global density function - for monitoring convergence
#'
#' @export
#' @method plot_global_density ndp
#' @param x ndp object
#' @param p probability mass contained in uncertainty interval
#' @return ggplot plot object
#'
plot_global_density.ndp <- function(x, p = 0.9, r = NULL,transform = TRUE ){

	if(p >= 1 || p <= 0 )
		stop("p must be in (0,1)")

    plt <- x$global_density %>%
		dplyr::mutate(Distance = (transform==TRUE)*Distance + (transform==FALSE)*qnorm(Distance / ceiling(max(x$if_df$Distance)))) %>%
		dplyr::group_by(Chain,Distance) %>%
        dplyr::summarise(lower = quantile(Global_Density,.5 - p / 2,na.rm=T),
                         med = median(Global_Density,na.rm=T),
                         upper = quantile(Global_Density,.5 + p /2,na.rm=T)) %>%
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
