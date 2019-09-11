
#' Random nonhomogenous poisson process generator
#' 
#' @param nsim number of processes to simulate
#' @param lambda intensity function (unnormalized)
#' @param max optional maximum value of intensity function
#' @param interval 1D interval overwhich process is realized
#' @param seed rng seed initializer
#' 
#' @export 
rnhpp <- function(nsim = 1, lambda = function(x) x^2, max = NULL, interval = c(0,1),
                  seed = NULL){
    if(!is.null(seed))
        set.seed(seed)
    if(is.null(max))
        stop("Must provide max")
    p_t <- function(x) lambda(x) / max
    out <- data.frame()
    i <- 1
    while(i <= nsim){
        ts <- numeric()
        t <- 0
        I <- 0
        S_I <- 0
        while(TRUE){
            u_1 <- runif(1)
            t <- t - log(u_1)/max
            if(t > interval[2])
                break
            u_2 <- runif(1)
            if(u_2 <= p_t(t) ){
                ts <- c(ts,t)
                I <- I + 1
                S_I <- t
            }
        }
        if(length(ts)!=0){
            out <- rbind(out,data.frame(sim_id = i,event_times = ts))
            i <- i + 1
        }
    }
    return(out)
}
