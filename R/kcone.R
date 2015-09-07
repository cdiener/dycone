# kcone.R Copyright 2015 Christian Diener <ch.diener@gmail.com> MIT license. See
# LICENSE for more information.

# Functions to analyze the k-cone

DC_SINGLECOL <- colorRampPalette(c("white", "darkred"))
DC_DIVCOL <- colorRampPalette(c("blue", "white", "tomato"))
DC_TRANSCOL <- colorRampPalette(c("seagreen", "darkgoldenrod1", "tomato"))

#' Calculates the euclidean norm of a vector
#'
#' @export
#' @param x A numeric vector.
#' @return The norm of the vector x.
#' @examples
#' x <- rnorm(10)
#' enorm(x)
enorm <- function(x) sqrt(sum(x * x))

#' Get the k-cone from a given flux cone.
#'
#' \code{kcone()} builds the k-cone from a given flux cone basis V and a set
#' of metabolic terms M by calculating \eqn{M^{-1}V}.
#' 
#' @seealso \code{\link{polytope_basis}} for a way to calculate the flux cone, 
#'  \code{\link{ma_terms}} to get mass action termsq 
#' @export
#' @keywords internal basis
#' @param V The flux cone basis. Rows denote the dimensions and columns the
#'  individual basis vectors.
#' @param m_terms The metabolic terms. A vector with as many entries as rows in V.
#'  All entries hsould be larger than 0.
#' @param normalize Whether the basis vectors should be scaled to unit length.
#'  Not recommened for differential analysis.
#' @return The stuff :O
#' @examples
#' data(eryth)
#' S <- matrix(c(1,0,0,1,-1,0, 0, -1), nrow=2)
#' rownames(S) <- c('A', 'B')
#' V <- polytope_basis(S)
#' K <- kcone(V, runif(ncol(S)))
kcone <- function(V, mats, normalize = FALSE) {
    K <- diag(1/mats) %*% V
    if (normalize) 
        K <- apply(K, 2, function(x) x/enorm(x))
    K <- as.matrix(K)
    class(K) <- append(c("basis", "kcone"), class(K))
    return(K)
}

#' Gets the null-space for the given k-cone assuming that all reactions
#' are irreversible.
#' 
#' @seealso \code{\link{polytope_basis}} for arbitrary reversibility 
#'  constraints.
#' @export
#' @keywords basis, k-cone
#' @param s_matrix The stochiometric matrix to be used.
#' @param m_terms The metabolic terms. 
#' @return A matrix containing the basis vectors in the columns.
#' @examples
#' S <- matrix(c(1,0,0,1,-1,0, 0, -1), nrow=2)
#' rownames(S) = c('A', 'B')
#' K <- kcone_null(S, runif(ncol(S)))
kcone_null <- function(s_matrix, m_terms) {
    SM <- s_matrix %*% diag(m_terms)
    
    basis <- MASS::Null(t(SM))
    class(basis) <- append("basis", class(basis))
    
    return(basis)
}

#' Gets the polytope basis of the k-cone.
#' 
#' @export
#' @keywords k-cone, basis
#' @param s_matrix The stochiometric matrix.
#' @param m_terms The metabolic terms for each reaction. Default to
#'  calculation of flux cone.
#' @return A matrix containing the basis vectors (normalized to a sum of 1) 
#'  in its columns.
#' @examples
#' S <- matrix(c(1,0,0,1,-1,0, 0, -1), nrow=2)
#' rownames(S) = c('A', 'B')
#' V <- polytope_basis(S)
polytope_basis <- function(s_matrix, m_terms = rep(1, ncol(s_matrix))) {
    mat <- m_terms
    const_matrix <- rcdd::d2q(-diag(ncol(s_matrix)))
    const_b <- rcdd::d2q(rep(0, ncol(s_matrix)))
    NC <- rcdd::d2q(as.matrix(s_matrix %*% diag(mat)))
    b <- rcdd::d2q(rep(0, nrow(s_matrix)))
    
    hp <- rcdd::makeH(const_matrix, const_b, NC, b)
    hp <- rcdd::addHeq(rep(1, ncol(s_matrix)), 1, hp)
    hp <- rcdd::redundant(hp)$output
    vrep <- rcdd::scdd(hp)
    vp <- vrep$output[, -(1:2)]
    
    basis <- as.matrix(apply(rcdd::q2d(t(vp)), 2, function(x) x/enorm(x)))
    class(basis) <- append("basis", class(basis))
    
    return(basis)
}

#' Translates the eigenvectors to the corresponding stability term.
#' 
#' @seealso \code{\link{stability_analysis}} for an analysis of the entire
#'  basis.
#' @export
#' @keywords stability
#' @param evs A vector of (complex) eigenvalues.
#' @return A single string of either stable, unstable, limit cycle or 
#'  zero jacobian.
#' @examples
#' print(as.stability(c(0+2i, 0-3i)))
as.stability <- function(evs) {
    evs[abs(evs) < sqrt(.Machine$double.eps)] <- 0
    
    r <- evs
    if (all(Re(r) == 0) && any(Im(r) != 0)) 
        return("limit cycle")
    if (all(Re(r) == 0)) 
        return("zero jacobian")
    if (all(Re(r) <= 0)) 
        return("stable") else return("unstable")
}

#' Performs a stability analysis for each basis vector in a k-cone.
#' 
#' @seealso \code{\link{polytope_basis}} and \code{}\link{kcone} to 
#'  calculate the k-cone.
#' @export
#' @keywords stability
#' @param basis A basis matrix, where each column denotes a basis vector.
#' @param s_matrix The stochiometric matrix of the model.
#' @param concs A named vector of metabolite concentrations.
#' @return A data frame containing the a string describing the type of 
#'  stability along with its eigenvalues for each vector in the basis.
#' @examples
#' data(eryth)
#' S <- matrix(c(1,0,0,1,-1,0, 0, -1), nrow=2)
#' rownames(S) = c('A', 'B')
#' V <- polytope_basis(S)
#' concs <- runif(nrow(S))
#' names(concs) <- rownames(S)
#' K <- kcone(V, ma_terms(S, concs))
#' stab <- stability_analysis(K, S, concs)
stability_analysis <- function(basis, s_matrix, concs) {
    J <- jacobian(s_matrix, concs)
    evs <- NULL
    stab <- NULL
    
    if (!is.null(ncol(basis))) {
        nc <- ncol(basis)
        data <- lapply(1:ncol(basis), function(i) {
            b <- basis[, i]/enorm(basis[, i])
            evs <- t(eigen(s_matrix %*% diag(b) %*% J)$values)
            stab <- as.stability(evs)
            return(data.frame(what = stab, evs = evs))
        })
    } else {
        nc <- 1
        b <- basis/enorm(basis)
        evs <- t(eigen(s_matrix %*% diag(b) %*% J)$values)
        stab <- as.stability(evs)
        data <- data.frame(what = stab, evs = evs)
    }
    
    if (nc > 1) 
        data <- do.call(rbind, data)
    res <- cbind(row.names = 1:nc, data)
    names(res)[1:nrow(s_matrix) + 1] <- paste0("ev", 1:nrow(s_matrix))
    
    return(res)
}


#' Calculates the occupation for a basis. Here, occupation denotes 
#' the faction of basis vectors for each varibale where the variable 
#' is not zero.
#' 
#' @seealso \code{\link{polytope_basis}} to get a k-cone.
#' @export
#' @keywords basis
#' @param basis A basis matrix, where each column denotes a basis vector.
#' @return A vector of length nrow(basis) containing the occupation for
#'  each dimension.
#' @examples
#' S <- matrix(c(1,0,0,1,-1,0, 0, -1), nrow=2)
#' rownames(S) = c('A', 'B')
#' V <- polytope_basis(S)
#' occ <- occupation(V)
occupation <- function(basis) {
    if (!("basis" %in% class(basis))) 
        stop("object is not a basis!")
    
    return(apply(basis, 1, function(r) sum(abs(r) > sqrt(.Machine$double.eps))/length(r)))
}

#' Plots a heatmap of the basis.
#'
#' @seealso \code{\link{plot_red}} for a geometric visualization.
#' @export
#' @keywords basis, plot
#' @param b The basis matrix, containing the basis vectors in the columns.
#' @param n_cl The number of clusters individual basis vectors are grouped 
#'  into. NA means no clustering is performed.
#' @param ... other arguments passed to pheatmap.
#' @examples
#' S <- matrix(c(1,0,0,1,-1,0, 0, -1), nrow=2)
#' rownames(S) <- c('A', 'B')
#' V <- polytope_basis(S)
#' plot_basis(V)
plot_basis <- function(b, n_cl = NA, ...) {
    clustered <- F
    
    if (!is.na(n_cl)) {
        write("Basis are very large. Reducing by clustering...", file = "")
        cl <- kmeans(t(b), centers = n_cl, iter.max = 100, nstart = 3)
        write(sprintf("Mean in-cluster distance: %g.", mean(sqrt(cl$withinss/cl$size))), 
            file = "")
        b <- t(cl$centers)
        clustered <- T
    }
    
    if (clustered) {
        ann <- data.frame(n = cl$size)
        pheatmap::pheatmap(b, border = NA, col = DC_SINGLECOL(101), labels_row = 1:nrow(b), 
            annotation_col = ann)
    } else {
        pheatmap::pheatmap(b, border = NA, col = DC_SINGLECOL(101), labels_row = 1:nrow(b), 
            ...)
    }
}

#' Plots a geometric visualization of the given k-cone reduced into
#' two dimensions.
#'
#' Plotting the reduced k-cone is a two step procedure. Where the k-cones
#' are first projected into the two most informative dimensions by PCA. In 
#' case the basis is very large this is followed by clustering of the basis 
#' vectors to reduce overlap in the projection. In the plot, arrows denote
#' the basis vectors, or 'extreme rays' of the respective k-cones and the 
#' shaded area denotes the interior of the k-cone into which all feasible
#' k vectors must fall.
#' 
#' @seealso \code{\link{plot_basis}} for heatmap of the basis.
#' @export
#' @keywords basis, plot
#' @param basis_list A list of basis to be reduced.
#' @param arrows Whether to draw the individual basis vectors or just the 
#'  interior of the cone.
#' @param col If NULL colors are automatically generated from a green to red
#'  palette. If not NULL must be of the same length as basis_list and denotes
#'  the colors used to draw the individual basis.
#' @param n_cl The numbers of clusters to be used. If there are more than 1000
#'  basis vectors per basis, clustering is always performed since calculation 
#'  of the convex hull (interior area of the cone) is extremely lengthy 
#'  otherwise.
#' @examples
#' S <- matrix(c(1,0,0,1,-1,0, 0, -1), nrow=2)
#' rownames(S) <- c('A', 'B')
#' V <- polytope_basis(S)
#' plot_red(list(V)) # A single basis is fine as long as it is a list
plot_red <- function(basis_list, arrows = TRUE, col = NULL, n_cl = NULL, r_names=NULL) {
    if (!("list" %in% class(basis_list))) 
        stop("basis_list must be a list!")
    
    b_rows <- sapply(basis_list, nrow)
    b_cols <- sapply(basis_list, ncol)
    
    if (length(basis_list) > 1 && sd(b_rows) > 0) 
        stop("All basis need to have the same dimension!")
    if (any(b_rows < 2)) 
        stop("Basis must be at least 2-dimensional!")
    n <- length(basis_list)
    
    if (!is.null(r_names) && (length(r_names) != max(b_rows) || 
		!is.character(r_names))) stop("r_names must be a character vector with as
		many entries as the maximum row number in basis_list!")
    
    all_basis <- do.call(cbind, basis_list)
    
    if (nrow(all_basis) == 2) {
        pca <- list(sdev = c(1, 1))
        red <- lapply(basis_list, t)
        if (!is.null(r_names)) load_names <- as.matrix(r_names)
    } else {
        pca <- prcomp(t(all_basis))
        red <- lapply(basis_list, function(b) predict(pca, t(b))[, 1:2])
        if (!is.null(r_names)) 
			load_names <- t(apply(pca$rotation[,1:2], 2, function(l) 
				r_names[order(-abs(l))[1:min(6,min(b_rows))]]))
    }
    
    if (max(b_cols) > 1000 || !is.null(n_cl)) {
        if (is.null(n_cl)) 
            n_cl <- floor(min(1000, sqrt(mean(b_cols)/2)))
        
        write(sprintf("Basis are very large. Reducing to %d clusters...", n_cl), 
            file = "")
        
        cl <- lapply(red, function(b) suppressWarnings(kmeans(b, centers = n_cl, 
            iter.max = 100, nstart = 3)))
        
        dists <- sapply(cl, function(x) mean(sqrt(x$withinss/x$size)))
        write(sprintf("Mean in-cluster distance: %g.", mean(dists)), file = "")
        
        red <- lapply(cl, function(x) x$centers)
    }
    
    energy <- pca$sdev
    write(sprintf("Information captured in projection: %f%%.", sum(energy[1:2])/sum(energy) * 
        100), file = "")
    
    rs <- apply(do.call(rbind, red), 2, range)
    
    if (!is.null(r_names)) {
		xl = paste0(load_names[1,], collapse=", ")
		yl = paste0(load_names[2,], collapse=", ")
	} else {
		xl = "PC 1"
		yl = "PC 2"
	}
    plot(NULL, xlim = rs[, 1], ylim = rs[, 2], xlab = xl, ylab = yl)
    
    if (is.null(col)) 
        pal <- DC_TRANSCOL(n) else if (length(basis_list) != length(col)) 
        stop("col mus be the same length as the basis list!") else pal <- col
    
    for (i in 1:n) {
        a <- apply(red[[i]], 1, angle)
        ro <- red[[i]][order(a), ]
        
        p <- rbind(ro, c(0, 0))
        red_ids <- rcdd::redundant(rcdd::makeV(rcdd::d2q(p)))$redundant
        polygon(p[-red_ids, ], border = NA, col = adjustcolor(pal[i], alpha.f = 0.2))
        if (arrows) {
            arrows(x0 = 0, y0 = 0, x1 = ro[, 1], y1 = ro[, 2], angle = 15, length = 0.05, 
                col = pal[i])
        }
    }
}

# Helper function to calculate angls between tow vectors
angle <- function(x, y = 0:1, clockwise = T) {
    theta <- acos(sum(x * y)/sqrt(sum(x * x) * sum(y * y)))
    if (clockwise && x[1] < 0) 
        theta <- 2 * pi - theta
    return(theta)
}

# Helper function to calculate the length ratio of two vectors
scaling <- function(x, y = 0:1) {
    s <- sum(x * x)/sum(y * y)
    return(sqrt(s))
}

# return( data.frame(id1=1:ncol(b1), id2=closest[1,], d=closest[2,]) ) }

#' Evaluates whether a given point lies within a given k-cone. Assumes a 
#' pointed cone (all reactions irreversible).
#' 
#' @export
#' @keywords basis
#' @param x A single point (vector) or multiple points (matrix with points as
#'  columns) to be checked.
#' @param s_matrix The stochiometric matrix of the k-cone to be used.
#' @param m_terms The metabolic terms of the k-cone.
#' @param tol The numerical accuracy of the check. Defaults to the square
#'  of the double accuracy.
#' @return A single boolean of vector of booleans indicating whether the
#'  point(s) lie within the k-cone.
#' @examples
#' S <- matrix(c(1,0,0,1,-1,0, 0, -1), nrow=2)
#' rownames(S) <- c('A', 'B')
#' 
#' # Check whether a random point falls within the flux cone
#' inside(runif(ncol(S)), S, rep(1,ncol(S))) # probably not true
inside <- function(x, s_matrix, m_terms, tol = sqrt(.Machine$double.eps)) {
    right <- s_matrix %*% diag(m_terms) %*% x
    
    if (is.null(dim(x))) 
        return(all(right > -tol) & all(abs(right) < tol))
    return(apply(right, 2, function(r) all(abs(r) > -tol)))
}

#' Gets the eigendynamics for a given k.cone basis.
#'
#' The eigendynamics for a given k-cone is the set of k vectors reproducing
#' the original k-cone with the most accuracy, using only k basis-vectors.
#' Here, only the first eigendynamics is guaranteed to be consistent with
#' the k-cone (all entries positive), all the additional eigendynamics can
#' be intepreted as fine tuning of the first eigendynamics.
#' 
#' @export
#' @keywords basis
#' @param basis The k-cone to be reduced.
#' @param n The number of eigendynamics to extract. Defaults to only the first.
#' @return A matrix containing the first k eigendynamics in its columns.
#' @examples
#' S <- matrix(c(1,0,0,1,-1,0, 0, -1), nrow=2)
#' rownames(S) <- c('A', 'B')
#' V <- polytope_basis(S)
#' ed <- eigendynamics(V) # gets the eigenpathways
eigendynamics <- function(basis, n = 1) {
    if (is.null(dim(basis))) 
        return(basis)
    
    eps <- svd(basis)$u
    eps[abs(eps) < .Machine$double.eps] <- 0
    if (any(eps[, 1] < 0)) 
        eps <- -eps
    
    return(eps[, 1:n])
}

#' Helper function to get the differential regulation between two approximations
#' of kinetic constansts.
#'
#' @seealso \code{\link{hyp}} for a complete differential analysis.
#' @export
#' @param k1 A vector containing kinetic constants >=0.
#' @param k2 A vector containing kinetic constants >=0.
#' @param reacts A reaction lists.
#' @return A data frame containig the reaction indices, reactions, kind of regulation
#'  and log2-fold changes of k2 relative to k1. Log-fold changes where either of the
#'  constants is 0 are evaluated to 0.
single_hyp <- function(k1, k2, reacts) {
    reacts <- make_irreversible(reacts)
    logfold <- log(k2, 2) - log(k1, 2)
    logfold[!is.finite(logfold)] <- 0
    
    
    idx <- id <- r_s <- kind <- fac <- NULL
    for (i in 1:length(logfold)) {
        r <- reacts[[i]]
        idx <- c(idx, i)
        id <- c(id, r$abbreviation)
        r_s <- c(r_s, to_string(r, name = F))
        kind <- c(kind, if (logfold[i] == 0) "same" else if (logfold[i] > 0) "up" else "down")
        fac <- c(fac, logfold[i])
    }
    
    res <- data.frame(idx = idx, name = id, reaction = r_s, type = kind, log2_fold = fac)
    
    return(res)
}

#' Identifies hypothesis for differentially regulated reactions between a set of
#' normal and disease conditions.
#'
#' \code{hyp()} is the main driver for differential analysis. It first uses a
#' reduction method to approximate the fold-change of the kinetic constants. This
#' is followed by calculation of all combinatorial log-fold changes within the 
#' normal group and between the disease and normal group. Finally, the wilcoxon
#' ranked sum statistics are obtained for each reaction and the p-values 
#' corrected.
#'
#' @export
#' @keywords hypothesis, k-cone, analysis
#' @param normal A matrix or data frame containing the metabolic terms of the 
#'  normal condition in the columns.
#' @param disease A matrix or data frame containing the metabolic terms of the 
#'  diseasse condition in the columns.
#' @param reacts The reaction list.
#' @param type The type of analysis to be performed. Either 'transformation',
#'  'optimization' or 'raw' for a pass-through option.
#' @param sorted Whether the results should be sorted by p-value and mean log-fold
#'  change.
#' @param v_opt Only needs to be set if type=='optimization'. Defines a weights
#'  specifying the optimization criterion. If F are the steady state fluxes of
#'  the specified system, the optimization criterion will be maximization of
#'  sum(v_opt*F). Must be the same length as there are irreversible reactions. 
#' @param correction_method A correction method for the multiple test p-values.
#'    Takes the same arguments as the method argument in p.adjust.
#' @param full If TRUE also returns the individual log fold changes along with
#'    the differential regulation data.
#' @return If full is FALSE only returns the generated hypothesis as a data 
#'  frame. If full is TRUE returns a list of generated hypothesis and the 
#'  individual log fold changes between all reference basis and between 
#'  reference and treatments. The full output will be a list with following 
#'  elements:
#'  \describe{
#'  \item{hyp}{The generated hypotheses together with statistics and reactions.}
#'  \item{lfc_normal}{The log2-fold changes of enzyme activity within the normal 
#'    group for each of the irreversible reactions. Those are never sorted, so
#'    the first entry corresponds to the first reaction, etc.}
#'  \item{lfc_disease}{The log2-fold changes of enzyme activity within between 
#'    the disease and normal group for each of the irreversible reactions. 
#'    Also never sorted.}
#'  \item{obj_normal}{Only if type=='optimization'. The value of the objective
#'    function for each of the metabolic terms in 'normal'.}
#'  \item{obj_disease}{Only if type=='optimization'. The value of the objective
#'    function for each of the metabolic terms in 'disease'.}
#'  }
hyp <- function(normal, disease, reacts, type = "transformation", correction_method = "fdr", 
    sorted = T, v_opt = NULL, full = FALSE) {
    # Create reference data
    cref <- combn(1:ncol(normal), 2)
    normal <- as.matrix(normal)
    disease <- as.matrix(disease)
    
    if (type == "transformation") {
        normal <- 1/normal
        disease <- 1/disease
    } else if (type == "optimization") {
        M <- cbind(normal, disease)
        S <- stochiometry(reacts)
        if (length(v_opt) != ncol(S) || !is.numeric(v_opt)) 
            stop("v_opt must be numeric and have the same length\n\t\t\tas there are irreversible reactions in reacts.")
        if (requireNamespace("foreach", quietly = TRUE)) {
            i <- 1:ncol(M)
            opt <- foreach::"%dopar%"(foreach::foreach(i = i, .combine = cbind), 
                dba(v_opt, S, M[, i], lower = 0, upper = 1))
        } else {
            opt <- lapply(1:ncol(M), function(i) dba(v_opt, S, M[, i], lower = 0, 
                upper = 1))
            opt <- do.call(cbind, opt)
        }
        normal <- opt[, 1:ncol(normal)]
        disease <- opt[, (ncol(normal) + 1):ncol(M)]
    } else if (type != "raw") 
        stop("type must be either 'transformation' or 'optimization' :(")
    
    res <- single_hyp(normal[, 1], normal[, 1], reacts)
    res <- res[, -ncol(res)]
    
    lfc_n <- apply(cref, 2, function(idx) {
        h <- single_hyp(normal[, idx[1]], normal[, idx[2]], reacts)
        return(h$log2_fold)
    })
    lfc_n <- cbind(lfc_n, -lfc_n)
    
    # Create differential analysis
    ctreat <- expand.grid(1:ncol(normal), 1:ncol(disease))
    
    lfc_d <- apply(ctreat, 1, function(idx) {
        h <- single_hyp(normal[, idx[1]], disease[, idx[2]], reacts)
        return(h$log2_fold)
    })
    
    # Generate statistics
    stats <- lapply(1:nrow(res), function(i) {
        n_data <- as.numeric(lfc_n[i, ])
        d_data <- as.numeric(lfc_d[i, ])
        test <- tryCatch({
            wilcox.test(x = d_data, conf.int = T)
        }, error = function(e) {
            list(p.value = 1, conf.int = rep(mean(n_data), 2))
        }, warning = function(w) {
            suppressWarnings(wilcox.test(x = d_data, y = n_data, conf.int = T))
        })
        
        return(data.frame(sd_normal = sd(n_data), sd_disease = sd(d_data), mean_log_fold = mean(d_data), 
            ci_low = test$conf.int[1], ci_high = test$conf.int[2], pval = test$p.value))
    })
    
    res <- cbind(res, do.call(rbind, stats))
    reg <- sapply(res$mean_log_fold, function(x) if (x == 0) 
        "same" else if (x > 0) 
        "up" else "down")
    res$type <- factor(reg)
    res$pval <- p.adjust(res$pval, method = correction_method)
    if (sorted) 
        res <- res[order(res$pval, -abs(res$mean_log_fold)), ]
    
    if (full) {
        lfc_n <- data.frame(normal = lfc_n)
        lfc_d <- data.frame(disease = lfc_d)
        res <- list(hyp = res, lfc_normal = lfc_n, lfc_disease = lfc_d)
        if (type == "optimization") {
            objval <- v_opt %*% (opt * M)
            res <- c(res, list(obj_normal = objval[1:ncol(normal)], 
				obj_disease = objval[-(1:ncol(normal))], k_normal = normal,
				k_disease = disease))
        }
    }
    
    return(res)
} 
