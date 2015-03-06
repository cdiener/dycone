#  approx.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.


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
approx_kinetics = function(reacts, concs, max_subs=1, degree=3, reversible=FALSE) {
	approx_ids = r_order(reacts) <= max_subs
	new_reacts = NULL
	new_concs = NULL
	
	for( i in 1:length(reacts) ) {
		if( approx_ids[i] ) {
			new_r = rep(reacts[[i]], degree+1)
			new_reacts = append( new_reacts, new_r )
			new_c = Ispline(concs[i], knots=c(0,concs[i],2*concs[i]), degree=degree)
			new_c = append(new_concs, new_c)
		}
		else
		{
			new_reacts = append( new_reacts, reacts[[i]] )
			new_concs = append( new_concs, concs[i] )
		}
	}
	
	S = get_stochiometry(new_reacts, reversible=reversible)
	
	res = list(S=S, concs=new_concs, degree=degree, approx=approx_ids)
	class(res) = append(class(res), "am")
}
