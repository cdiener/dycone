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
dba = function(a, s_matrix, v_terms, lower=0, upper=20) {
	mat = v_terms
	const_matrix = rcdd::d2q( rbind(-diag(ncol(s_matrix)), diag(ncol(s_matrix))) )
	const_b = rcdd::d2q(rep(c(-lower,upper),each=ncol(s_matrix)))
	NC = rcdd::d2q( as.matrix( s_matrix%*%diag(mat) ) )
	b = rcdd::d2q(rep(0,nrow(s_matrix)))
	
	hp = rcdd::makeH(const_matrix, const_b, NC, b)
	
	sol = rcdd::lpcdd(hp, rcdd::d2q(a), minimize=FALSE)
	if(sol$solution.type!="Optimal") stop("Optimization is inconsistent!")
	
	return(rcdd::q2d(sol$primal.solution))
}

#' Finds all reaction indices using the given subtrates and products.
#'
#' @param r A reaction list.
#' @param S a list of substrates to be searched for or NULL.
#' @param P a list of products to be searched for or NULL.
which_reaction = function(r, S, P) {
	have_it = sapply(r, function(x) any(S %in% x$S) || any(P %in% x$P))
	
	return(which(have_it))
}
