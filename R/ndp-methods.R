#' Methods for ndp objects
#' 
#' Common methods that allow for easy model interrogation
#' 
#' @name ndp-methods
#' @aliases nsamples 
#' 
#' @param object ndp object
#' @param ... Ignored
#' 

#'
#' @rdname ndp-methods
#' @export
#' @importFrom rstantools nsamples
nsamples.ndp <- function(object, ...) {
  nrow(object$alpha)
}




# Internal ------------------------------------

