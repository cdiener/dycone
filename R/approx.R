#  approx.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.


split_react = function(r, n)
{
	res = list()
	for(i in 1:n) {
		new_r = r
		new_r$S = sapply(new_r$S, function(x) paste0(x, "_", i))
		names(new_r$S) = NULL
		res = append(res, list(new_r))
	}
	
	return(res)
}

#' Approximates the kinetics of low order reactions with splines allowing the
#' mimicry of arbitrary non-linear kinetics.
#'
#' @param reacts The reactions to be approximated
#' @param concs Vector of steady state concentrations
#' @param max_subs Only reactions with at most max_subs different substrates will be
#' 		approximated.
#' @param degree Polynomial degree of the approximation
#' @param reversible Whether reversible reactions are allowed or will be split.
#' @return An R object conatining the following elements:
#' 	\describe{
#'		\item{S}{The stochiometric matrix of the approximated system. Here
#'			the matrix will contain spoof reactions in order to realize the
#'			representation}
#'		\item{concs}{The approximated concentrations. For approximated reactions those
#'			are the values of the spline basis function at steady state concenyration.}
#'		\item{degree}{The polynomial degree of the approximation.}
#'		\item{approx}{Boolean vector indicating for each of the original reactions whether
#'				it has been approximated.}
#'}
approx_kinetics = function(reacts, concs, max_subs=1, degree=3, 
					reversible=FALSE, const="none") {
	approx_ids = r_order(reacts) <= max_subs & !constant_flux(reacts)
	new_reacts = NULL
	concs_idx = rep(FALSE, length(concs))
	names(concs_idx) = names(concs)
	new_concs = NULL
	
	for( i in 1:length(reacts) ) {
		if( approx_ids[i] ) {
			new_r = split_react(reacts[[i]], degree+1)
			new_reacts = append( new_reacts, new_r )
			concs_idx[reacts[[i]]$S] = TRUE
		}
		else
		{
			new_reacts = append( new_reacts, list(reacts[[i]]) )
		}
	}
	class(new_reacts) = append(class(new_reacts), "reactions")
	
	for( i in 1:length(concs) ) {
		if( concs_idx[i] ) {
			new_c = Ispline(concs[i], degree, c(0,concs[i],2*concs[i]) )
			names(new_c) = paste0(names(concs)[i],"_",1:(degree+1))
			new_concs = c(new_concs, concs[i], new_c)
		}
		else new_concs = c(new_concs, concs[i])
	}
	print(new_reacts)
	
	S = get_stochiometry(new_reacts, reversible=reversible, const=const)
	
	res = list(S=S, concs=new_concs, degree=degree, r_idx=approx_ids, s_idx=concs_idx)
	class(res) = append(class(res), "am")
	
	return(res)
}
