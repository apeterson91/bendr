#' Retrieve Grouped Data Structure
#'
#' @export
#' @param x benvo object
#' @param formula bendr formula
setGeneric("groupvo",function(x,formula) standardGeneric("groupvo")) 


#' Grouped
#' 
#' @describeIn groupvo return list of grouped data structures 
#' @export
#'
setMethod("groupvo","Benvo",function(x,formula)){

	id <- formula.tools::lhs(formula)
	vars <- formula.tools::rhs.vars(formula)
	if(length(vars)>1){
		stop("This function only takes 1 BEF as an argument")
	}

	R <- max(x@bef_data[[x]]$Distance)

	r <- x@bef_data[[x]] %>% 
		dplyr::arrange(!!id) %>% 
		dplyr::select(Distance) %>%  ## return scaled distances
		dplyr::mutate(Distance = qnorm(Distance/max(Distance))) %>%
			dplyr::pull()

	n <- x@bef_data[[x]] %>% 
		dplyr::arrange({{id}}) %>%
		dplyr::group_by({{dplyr::id_col}}) %>%
		dplyr::count() %>% 
		dplyr::ungroup() %>% dplyr::mutate(start = cumsum(n)) %>% 
		dplyr::mutate(start_ = tidyr::replace_na(dplyr::lag(start),0)) %>% 
		dplyr::select(-start) %>% dplyr::rename(start = start_,go=N) %>% 
		dplyr::select(start,go) %>% as.matrix()
	J <- nrow(n)

	return(list(r = r,n = n, J = J, R = R))

}


# for internal use ----------------------------------------------------------------------------





