#  linprog.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

#' Maximizes a given flux within the k-cone.
#'
#' @param a A vector a such that the fluxes a'%*%v are maximized.
#' @param s_matrix The stochiometrix matrix to be used (must be irreversible).
#' @param v_terms The corresponding flux products.  
dba = function(a, s_matrix, v_terms) {
	mat = v_terms
	const_matrix = rcdd::d2q( -diag(ncol(s_matrix)) )
	const_b = rcdd::d2q(rep(0,ncol(s_matrix)))
	NC = rcdd::d2q( as.matrix( s_matrix%*%diag(mat) ) )
	b = rcdd::d2q(rep(0,nrow(s_matrix)))
	
	hp = rcdd::makeH(const_matrix, const_b, NC, b)
	hp = rcdd::addHeq(rep(1,ncol(s_matrix)), 1.0, hp)
	
	sol = rcdd::lpcdd(hp, rcdd::d2q(a), minimize=FALSE)
	print(sol)
	if(sol$solution.type!="Optimal") stop("Optimization is inconsistent!")
	
	return(sol)
}
