#' Plot method for ndp objects
#'
#' The \code{plot} method for \link{ndp-objects} provides convenient 
#' access to several types of plots useful for visualizing
#' ndp models
#'
#' @method plot ndp 
#' @export
#' @param  x an ndp object
#' @importFrom ggplot2 ggplot aes_string xlab %+replace% theme
#'
#'
plot.ndp <- function(x, plotfun = "intervals", pars = NULL,
                         regex_pars = NULL, ...) {
  
    stop("Not yet implemented")
}

#' @export
plot_pairs <- function(x)
    UseMethod("plot_pairs")

#' @export
plot_intensities <- function(x, p = .9, pi_threshold = .1)
    UseMethod("plot_intensities")

#' @export
plot_traceplots <- function(x, par="alpha")
    UseMethod("plot_traceplots")

#' plots pairwise probability clustering plot
#' @method plot_pairs ndp
#' @export 
#' @param x ndp object
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn
#' 
plot_pairs.ndp <- function(x){

    p <- dplyr::as_tibble(x$pmat) %>% dplyr::mutate(Group_1 = 1:dplyr::n()) %>%
        tidyr::gather(dplyr::contains("V"), key = "Group_2", value = "Probability") %>% 
        dplyr::mutate("Group_2" = as.numeric(stringr::str_replace(Group_2,"V",""))) %>%
        dplyr::filter(Group_1 > Group_2) %>%
        ggplot(aes(x=Group_1,y=Group_2,fill=Probability)) + 
        geom_tile() + scale_fill_gradientn(colours=rainbow(10),limits=c(0,1)) + ggplot2::theme_bw() + ggplot2::labs(title = "Pairwise Probability of Function Clustering")

    return(p)
}

#' Plots intensity function densities
#' @export
#' @method plot_intensities ndp
#' @param x ndp object
#' @param p the probability for the credible interval
#' @param pi_threshold intensities with probability of assignment greater than pi_threshold are plotted 
#' default is zero
#' 
plot_intensities.ndp <- function(x, p = .9,pi_threshold = .1){

    ks_to_keep <- which(apply(as.matrix(x$pi),2,median)>pi_threshold)

    p <- x$if_df %>% dplyr::group_by(Intensity_Function,Distance) %>%
        dplyr::summarise(lower = quantile(Density,.5 + p/2,na.rm=T),
                         med = median(Density,na.rm=T),
                         upper = quantile(Density,.5 + p /2,na.rm=T)) %>%
    dplyr::filter(Intensity_Function %in% ks_to_keep) %>%
    ggplot(aes(x=Distance,y=med)) + ggplot2::geom_line() + 
    ggplot2::geom_ribbon(aes(ymin= lower,ymax=upper),alpha=0.3) + 
    ggplot2::theme_bw() + 
    ggplot2::facet_wrap(~Intensity_Function) + 
    ggplot2::theme(strip.background = ggplot2::element_blank()) + 
    ggplot2::labs(title = "Intensity Functions",
         subtitle = paste0("Shaded area indicates ",p,"% Credible Interval"),
         y = "Density")
    return(p)
}


#' Plots traceplots of various groups of parameters
#' @export
#' @method plot_traceplots ndp
#' @param x ndp object
#' 
plot_traceplots.ndp <- function(x,par="alpha"){

    p <- ggmcmc::ggs_traceplot(ggmcmc::ggs(x[[par]]))

    return(p)
}
