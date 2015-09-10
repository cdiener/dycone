# linprog.R Copyright 2015 Christian Diener <ch.diener@gmail.com> MIT license.
# See LICENSE for more information.

# Helper function to convert various input formats for the objective reaction
# Also performs a set of validity checks
build_objective <- function(o, S) {
    if (is.list(o)) {
        if (!all(c("S", "P", "N_S", "N_P") %in% names(o)))
            stop("objective list must contain entries for S, P, N_S and N_P.")
        l <- sapply(c("S", "P", "N_S", "N_P"), function(i) length(o[[i]]))
        if (l[1] - l[3] != 0) 
            stop("S and N_S must have the same length")
        if (l[2] - l[4] != 0) 
            stop("P and N_P must have the same length")
        valid_S <- !is.na(o$S)
        valid_P <- !is.na(o$P)
        sto <- c(-o$N_S[valid_S], o$N_P[valid_P])
        mets <- c(o$S[valid_S], o$P[valid_P])
        col <- rep(0, nrow(S))
        names(col) <- rownames(S)
        col[mets] <- sto
    } else if (is.vector(o)) {
        if (is.null(names(o))) {
            if (length(o) != nrow(S)) 
                stop("Unnamed objective vector must have the same length as rows in S.")
            col <- o
        } else {
            col <- rep(0, nrow(S))
            names(col) <- rownames(S)
            col[names(o)] <- o
        }
    } else stop("objective must be a list or vector.")
    
    return(col)
}

# Helper function to allow single or multi-value constraints
expand_constraints = function(v, n) {
    if (length(v) == n) return(v)
    else if (length(v) == 1) return(rep(v,n))
    else stop("There must be a single flux constraint or one for each reaction.")
}

#' Maximizes the flux through a given objective reaction.
#'
#' @param obj The objective reaction, whose flux will be be maximized. Can 
#'  be any of the following three:
#'  \itemize{
#'  \item{A list containing at least the four vectors S, P, N_S, and N_P, which
#'      which contain the names of substrates and products and their respective
#'      stochiometries.}
#'  \item{A unnamed vector containing the corresponding column in the 
#'      stochiometric matrix.}
#'  \item{A named vector containing only the non-zero entries in the respective
#'      column of the stochiometric matrix bein named by their respective 
#'      metabolite.}
#'  }
#' @param S The stochiometrix matrix to be used (must be irreversible).
#' @param v_min Lower bounds for the reaction fluxes. Can be a single value or a 
#'  vector containing one value for each reaction.
#' @param v_max Upper bounds for the reaction fluxes. Can be a single value or a 
#'  vector containing one value for each reaction.
#' @return A vector with the same length as columns in \code{S} containing the 
#'  solution of the optimization.
#' @export
fba <- function(obj, S, v_min = 0, v_max = 1) {
    v_min = expand_constraints(v_min, ncol(S))
    v_max = expand_constraints(v_max, ncol(S))
    
    S <- cbind(S, build_objective(obj, S))
    a <- rep(0, ncol(S))
    a[ncol(S)] <- 1
    
    const_matrix <- rcdd::d2q(rbind(-diag(ncol(S)), diag(ncol(S))[-ncol(S),]))
    const_b <- rcdd::d2q(c(-v_min, 0, v_max))
    NC <- rcdd::d2q(as.matrix(S))
    b <- rcdd::d2q(rep(0, nrow(S)))
    
    hp <- rcdd::makeH(const_matrix, const_b, NC, b)
    
    sol <- rcdd::lpcdd(hp, rcdd::d2q(a), minimize = FALSE)
    if (sol$solution.type != "Optimal") 
        stop(sprintf("Error in optimization! (%s)", sol$solution.type))
    
    return(rcdd::q2d(sol$primal.solution))
}

#' Maximizes a set of given fluxes within the k-cone using parsimonous FBA.
#' Here, 'parsimonous' means that from all the possible maxima we will return
#' the one minimizing the overall flux sum, thus, being the most 'economic' one
#' for the organism.
#'
#' @param obj The objective reaction, whose flux will be be maximized. Can 
#'  be any of the following three:
#'  \itemize{
#'  \item{A list containing at least the four vectors S, P, N_S, and N_P, which
#'      which contain the names of substrates and products and their respective
#'      stochiometries.}
#'  \item{A unnamed vector containing the corresponding column in the 
#'      stochiometric matrix.}
#'  \item{A named vector containing only the non-zero entries in the respective
#'      column of the stochiometric matrix bein named by their respective 
#'      metabolite.}
#'  }
#' @param S The stochiometrix matrix to be used (must be irreversible).
#' @param v_min Lower bounds for the reaction fluxes. Can be a single value or a 
#'  vector containing one value for each reaction.
#' @param v_max Upper bounds for the reaction fluxes. Can be a single value or a 
#'  vector containing one value for each reaction.
#' @return A vector with the same length as columns in \code{S} containing the 
#'  solution of the optimization.
#' @export
pfba <- function(obj, S, v_min = 0, v_max = 1) {
    v_min = expand_constraints(v_min, ncol(S))
    v_max = expand_constraints(v_max, ncol(S))
    
    S <- cbind(S, build_objective(obj, S))
    a <- rep(0, ncol(S))
    a[ncol(S)] <- 1
    
    const_matrix <- rcdd::d2q(rbind(-diag(ncol(S)), diag(ncol(S))[-ncol(S),]))
    const_b <- rcdd::d2q(c(-v_min, 0, v_max))
    NC <- rcdd::d2q(as.matrix(S))
    b <- rcdd::d2q(rep(0, nrow(S)))
    
    hp <- rcdd::makeH(const_matrix, const_b, NC, b)
    
    sol <- rcdd::lpcdd(hp, rcdd::d2q(a), minimize = FALSE)
    if (sol$solution.type != "Optimal") 
        stop(sprintf("Error in initial optimization! (%s)", sol$solution.type))
    
    keep_idx <- which(a != 0)    
    keep_val <- sol$primal.solution[keep_idx]
    a <- diag(ncol(S))[keep_idx,]
    pp <- rcdd::addHeq(a, keep_val, hp)
    sol <- rcdd::lpcdd(pp, rcdd::d2q(rep(1,ncol(S))), minimize = TRUE)
    
    if (sol$solution.type != "Optimal") 
        stop(sprintf("Error in parsimonous optimization! (%s)", sol$solution.type))
    
    return(rcdd::q2d(sol$primal.solution))
}

#' Finds all reaction indices using the given subtrates and products.
#'
#' @export
#' @param r A reaction list.
#' @param S a list of substrates to be searched or empty string.
#' @param P a list of products to be searched or empty string.
#' @return A data frame with two columns, the first specifying the index
#'  of the reaction and the second the metabolites found for that reaction.
which_reaction <- function(r, S = "", P = "") {
    have_it <- sapply(r, function(x) any(S %in% x$S) || any(P %in% x$P))
    idx <- which(have_it)
    out <- lapply(idx, function(i) {
        subs <- S[S %in% r[[i]]$S]
        prods <- P[P %in% r[[i]]$P]
        data.frame(idx = i, metabolites = c(subs, prods))
    })
    
    return(do.call(rbind, out))
}

#' Finds the closest point in within the k-cone to a query point p.
#' 
#' @export
#' @seealso \code{\link{fba}} for optimization within the k-cone.
#' @param p The query point. Can be outside or inside the k-cone.
#' @param S The stochiometric matrix.
#' @param m_terms The metbaolic terms to be used. Must be have 
#'  \code{ncol(S)} elements.
#' @return A vector containing the closest point to p within the k-cone.
closest <- function(p, S, m_terms) {
    if (!requireNamespace("quadprog", quietly = TRUE)) 
        stop("This function requires the quadprog package.")
    
    if (length(p) != ncol(S)) 
        stop("p does not have the correct dimension!")
    dp <- enorm(p)
    
    A <- S %*% diag(m_terms)
    A <- t(rbind(A, diag(ncol(S))))
    qp <- quadprog::solve.QP(diag(ncol(S)), p/dp, A, meq = nrow(S))
    
    
    return(qp$solution * dp)
} 
