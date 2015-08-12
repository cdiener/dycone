#  linprog.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

#' Maximizes a given flux within the k-cone.
#'
#' @param a A vector a such that the fluxes \eqn{a^Tv} are maximized.
#' @param s_matrix The stochiometrix matrix to be used (must be irreversible).
#' @param m_terms The corresponding flux products.
#' @param lower Lower bound for the reaction constants.
#' @param upper Upper bound for the reactions constants.
#' @return A vector with the same length as \code{a} containing the solution
#'  of the optimization.
#' @export
dba = function(a, s_matrix, m_terms, lower=0, upper=1) {
	const_matrix = rcdd::d2q( rbind(-diag(ncol(s_matrix)), diag(ncol(s_matrix))) )
	const_b = rcdd::d2q(c(-lower*m_terms,upper*m_terms))
	NC = rcdd::d2q( as.matrix( s_matrix%*%diag(m_terms) ) )
	b = rcdd::d2q(rep(0,nrow(s_matrix)))
	
	hp = rcdd::makeH(const_matrix, const_b, NC, b)
	
	sol = rcdd::lpcdd(hp, rcdd::d2q(a), minimize=FALSE)
	if(sol$solution.type!="Optimal") stop("Optimization is inconsistent!")
	
	return(rcdd::q2d(rcdd::qdq(sol$primal.solution, rcdd::d2q(m_terms))))
}

#' Finds all reaction indices using the given subtrates and products.
#'
#' @param r A reaction list.
#' @param S a list of substrates to be searched or empty string.
#' @param P a list of products to be searched or empty string.
#' @export
which_reaction = function(r, S="", P="") {
	have_it = sapply(r, function(x) any(S %in% x$S) || any(P %in% x$P))
	idx = which(have_it)
    out = lapply(idx, function(i) {
        subs = S[S %in% r[[i]]$S]
        prods = P[P %in% r[[i]]$P]
        data.frame(idx=i, metabolites=c(subs,prods)) })
    
	return(do.call(rbind, out))
}

closest = function(p, S, m_terms) {
    if (!requireNamespace("quadprog", quietly = TRUE))
        stop("This function requires the quadprog package.")
    
    if(length(p) != ncol(S)) stop("p does not have the correct dimension!")
    dp = enorm(p)
    
    A = S%*%diag(m_terms)
    A = t(rbind(A, diag(ncol(S))))
    qp = quadprog::solve.QP(diag(ncol(S)), p/dp, A, meq=nrow(S))

	
	return(qp$solution*dp)
}
