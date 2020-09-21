#' Retrieve Grouped Data Structure
#'
#' @export
#' @param x benvo object
#' @param formula bendr formula
groupvo <- function(x,formula) UseMethod("groupvo")


#' @describeIn groupvo return list of grouped data structures 
#' @export
groupvo.benvo <- function(x,formula){


	start_ <- start <- Distance <- go <- NULL
	id <- formula.tools::lhs(formula)
	vars <- formula.tools::rhs.vars(formula)

	if(length(vars)>1)
		stop("This function only takes 1 BEF as an argument")

	ix <- which(rbenvo::bef_names(x) == vars)
	
	if(!length(ix)){
		st <- glue::glue("{vars} is not included as a member of this Benvo")
		stop(st)
	}

	R <- ceiling(max(x$bef_data[[ix]]$Distance))

	r <- x$bef_data[[ix]] %>% 
		dplyr::arrange({{id}}) %>% 
		dplyr::select(Distance) %>%  ## return scaled distances
			dplyr::pull()

	n <- x$bef_data[[ix]] %>% 
		dplyr::arrange({{id}}) %>%
		dplyr::group_by({{id}}) %>%
		dplyr::count() %>% 
		dplyr::ungroup() %>% dplyr::mutate(start = cumsum(n)) %>% 
		dplyr::mutate(start_ = tidyr::replace_na(dplyr::lag(start),0)) %>% 
		dplyr::select(-start) %>% dplyr::rename(start = start_,go=n) %>% 
		dplyr::select(start,go) %>% as.matrix()
	J <- nrow(n)

	return(list(r = r,n = n, J = J, R = R))

}

