#  kcone.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

# Functions to analyze the k-cone

DC_SINGLECOL = colorRampPalette(c("white","darkred"))
DC_DIVCOL = colorRampPalette(c("blue", "white", "tomato"))
DC_TRANSCOL = colorRampPalette(c("seagreen", "darkgoldenrod1", "tomato"))

norm = function(x) sqrt(sum(x*x))

#' Get the k-cone from a given flux cone.
#'
#' \code{kcone()} builds the k-cone from a given flux cone basis V and a set
#' of metabolic terms M by calculating \eqn{M^{-1}V}.
#' 
#' @seealso \code{\link{polytope}} for a way to calculate the flux cone, 
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
#' V = polytope(eryth)
#' K = kcone(V, mats)
kcone = function(V, mats, normalize=FALSE) {
    K = diag(1/mats)%*%V
	if(normalize) K = apply(K,2,function(x) x/norm(x))
    K = as.matrix(K)
    class(K) = append(c("basis", "kcone"), class(K))
	return(K)
}

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
get_kcone_basis = function(s_matrix, v_terms) {	
	mat = v_terms
	SM = s_matrix %*% diag(mat)
	
	basis = MASS::Null( t(SM) )
	class(basis) = append("basis", class(basis))
	
	return( basis )
}

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
get_polytope_basis = function(s_matrix, v_terms=rep(1, ncol(s_matrix))) {
	mat = v_terms
	const_matrix = rcdd::d2q(-diag(ncol(s_matrix)))
	const_b = rcdd::d2q(rep(0,ncol(s_matrix)))
	NC = rcdd::d2q( as.matrix( s_matrix%*%diag(mat)) )
	b = rcdd::d2q(rep(0,nrow(s_matrix)))
	
	hp = rcdd::makeH(const_matrix, const_b, NC, b)
	hp = rcdd::addHeq(rep(1,ncol(s_matrix)), 1.0, hp)
	hp = rcdd::redundant(hp)$output
	vrep = rcdd::scdd(hp)
	vp = vrep$output[,-(1:2)]
	
	basis = as.matrix(apply(rcdd::q2d(t(vp)), 2, function(x) x/norm(x)))
	class(basis) = append("basis", class(basis))
	
	return( basis)
} 

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
get_stability = function(evs) {	
	evs[ abs(evs)<sqrt(.Machine$double.eps) ] = 0
	
	r = evs
	if ( all(Re(r)==0) && any(Im(r)!=0) ) return( "limit cycle" )
	if ( all(Re(r)==0) ) return ("zero jacobian")
	if ( all(Re(r)<=0) ) return( "stable" )
	else return ("unstable")
}

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
stability_analysis = function(basis, s_matrix, concs) {
	J = get_jacobian(s_matrix, concs)
	evs = NULL
	stab = NULL
	
	if(!is.null(ncol(basis))) {
		nc = ncol(basis)
		for( i in 1:ncol(basis) )
		{
			ev = eigen(s_matrix %*% diag(basis[,i]) %*% J)$values
			evs = rbind(evs, ev)
			stab = c(stab, get_stability(ev))
		}
	}
	else {
		nc = 1
		evs = t(eigen(s_matrix %*% diag(basis) %*% J)$values)
		stab = get_stability(evs)
	}
		
	res = data.frame(what=stab, ev=evs, row.names=1:nc)
	names(res)[1:nrow(s_matrix)+1] = paste0("ev",1:nrow(s_matrix))
	
	return( res )
}


#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
occupation = function(basis) {
	if( !("basis" %in% class(basis)) ) stop("object is not a basis!")
	
	return(apply(basis, 1, function(r) 
		sum(r>sqrt(.Machine$double.eps))/length(r)))
}

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
plot.basis = function(b, n_cl=sqrt(ncol(b)/2), ...) {	
    clustered = F
    
    if((ncol(b)>1e3 || n_cl != sqrt(ncol(b)/2)) && !is.na(n_cl)) {
        write("Basis are very large. Reducing by clustering...", file="")
        cl = kmeans(t(b), centers=n_cl, iter.max=100, nstart=3)
        write(sprintf("Mean in-cluster distance: %g.",
            mean(sqrt(cl$withinss/cl$size))), file="")
        b = t(cl$centers)
        clustered = T
    }
    
    if(clustered) {
        ann = data.frame(n=cl$size)
        pheatmap::pheatmap(b, border=NA, col=DC_SINGLECOL(101), 
            labels_row=1:nrow(b), annotation_col=ann)
    }
    else {
        pheatmap::pheatmap(b, border=NA, col=DC_SINGLECOL(101), 
            labels_row=1:nrow(b), ...)
    }
}

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
plot_red = function(basis_list, arrows=TRUE, col=NULL, n_cl=NULL) {
	if( !("list" %in% class(basis_list)) ) stop("basis_list must be a list!")
	
	b_rows = sapply(basis_list, nrow)
    b_cols = sapply(basis_list, ncol)
    
    if(length(basis_list)>1 && sd(b_rows)>0) 
        stop("All basis need to have the same dimension!")
    if(any(b_rows<2)) stop("Basis must be at least 2-dimensional!")
    n = length(basis_list)

    all_basis = do.call(cbind, basis_list)
    
    if(nrow(all_basis)==2) {
        pca = list(sdev=c(1,1))
        red = lapply(basis_list, t)
    }
    else {
        pca = prcomp(t(all_basis))
        red = lapply(basis_list, function(b) predict(pca, t(b))[,1:2])
    }
    
    if(max(b_cols)>1e3 || !is.null(n_cl)) {
        write("Basis are very large. Reducing by clustering...", file="")
        if(is.null(n_cl)) n_cl = sqrt(mean(b_cols)/2)
        cl = lapply(red, function(b) 
            kmeans(b,centers=n_cl, iter.max=100, nstart=3))
        
        dists = sapply(cl, function(x) mean(sqrt(x$withinss/x$size)))
        write(sprintf("Mean in-cluster distance: %g.",mean(dists)), file="")
        
        red = lapply(cl, function(x) x$centers)
    }
    
	rs = apply(do.call(rbind,red), 2, range)
	plot(NULL, xlim=rs[,1], ylim=rs[,2], xlab="PC 1", ylab="PC 2")
	
	if(is.null(col)) pal = DC_TRANSCOL(n)
	else if(length(basis_list)!=length(col)) 
		stop("col mus be the same length as the basis list!")
	else pal = col
	
	for( i in 1:n) {
		a = apply(red[[i]], 1, angle)
		ro = red[[i]][order(a),]
		
		p = rbind(ro,c(0,0))
		red_ids = rcdd::redundant( rcdd::makeV(rcdd::d2q(p)) )$redundant
		polygon(p[-red_ids,], border=NA, col=adjustcolor(pal[i], alpha.f=0.2))
		if(arrows) {
			arrows(x0=0, y0=0, x1=ro[,1], y1=ro[,2], angle=15, 
				length=0.05, col=pal[i])
		}
	}
	energy = pca$sdev^2
	
	write(sprintf("Total energy explained: %f%%.",
				sum(energy[1:2])/sum(energy)*100), file="")
}

angle = function(x, y=0:1, clockwise=T) {
	theta = acos( sum(x * y) / sqrt(sum(x * x) * sum(y * y)) ) 
	if(clockwise && x[1]<0) theta = 2*pi - theta 
	return(theta)
}

scaling = function(x, y=0:1) {
    s = sum(x * x)/sum(y * y)
    return(sqrt(s))
}

#basis_map = function(b1, b2) {
#	n1 = ncol(b1)
#	n2 = ncol(b2)
#	progress = F
#	cur = 0
	
#	if(n1*n2>1e5) {
#		write("Basis are pretty large. This might take a long time!", file="")
#		progress = T
#	}
#	closest = apply(b1, 2, function(x) {
#							d = apply(b2, 2, function(y) dist(rbind(x,y)))
#							i = which.min(d)
#							if(progress) {
#								cur <<- cur+1
#								cat("\r") 
#								cat(sprintf("Finished %.2f%%...",cur/n1*100))
#							}
#							return( c(i,d[i]) )
							
#							 })
#	cat("\n")
	
#	return( data.frame(id1=1:ncol(b1), id2=closest[1,], d=closest[2,]) )
#}

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
inside = function(x, s_matrix, v_terms, tol=sqrt(.Machine$double.eps)) {
	right = s_matrix%*%diag(v_terms)%*%x
	
	if(is.null(dim(x))) return(all(right>-tol) & all(abs(right)<tol))
	return( apply(right, 2, function(r) all(r>=0) & all(abs(r)<tol)) )
}

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
eigendynamics = function(basis, n=1) {
	if(is.null(dim(basis))) return(basis)
	
	eps = svd(basis)$u
	eps[ abs(eps)<.Machine$double.eps ] = 0
	if( any(eps[,1]<0) ) eps = -eps
	
	eps = apply(eps,2,function(x) x/sum(x))
	
	return(eps[,1:n])
}

#r_means = function(x) {
#	if(!is.null(dim(x))) x = rowMeans(x)
	
#	return(x)
#}

single_hyp = function(b1, b2, reacts, log_fold_tol=0) {
	reacts = make_irreversible(reacts)
	logfold = log(b2,2) - log(b1,2)
	logfold[!is.finite(logfold)] = 0
	
	
	idx = id = r_s = kind = fac = NULL
	for( i in 1:length(logfold)) {
		r = reacts[[i]]
		idx = c(idx, i)
		id = c(id, r$abbreviation)
		r_s = c(r_s, to_string(r, name=F))
		kind = c(kind, if(logfold[i]==0) "same" else if(logfold[i]>0) "up" else "down")
		fac = c(fac,logfold[i])
	}
	
	res = data.frame(idx=idx, name=id, reaction=r_s, type=kind, log2_fold=fac)
	
	return(res)
}

mabs_diff = function(x,y) {
	d = sapply(x, function(xi) sapply(y, function(yi) xi-yi))
	m = mean(abs(d), na.rm=T)
	if(is.na(m)) m = 0
	
	return(m)
}

#' Identifies hypothesis for differentially regulated reactions between a set of
#' normal and disease conditions. 
#'
#' @param normal A matrix or data frame containing the metabolic terms of the 
#'  normal condition in the columns.
#' @param disease A matrix or data frame containing the metabolic terms of the 
#'  diseasse condition in the columns.
#' @param type The type of analysis to be performed. Either 'transformation',
#'  'optimization' or 'raw' for a pass-through option.
#' @param sort Whether the results should be sorted by p-value and mean log-fold
#'  change.
#' @param v_opt A vector defining the optimization criterion. Must have the same
#'  length as there are irreversible reactions. 
#' @param correction_method A correction method for the multiple test p-values.
#'	Takes the same arguments as the method argument in p.adjust.
#' @param full If TRUE also returns the individual log fold changes along with
#'	the differential regulation data.
#' @return If full is TRUE returns a list of generated hypothesis and the individual
#' 	log fold changes between all reference basis and between reference and treatments.
#'	If full is FALSE only returns the generated hypothesis. 
#' @export
hyp = function(normal, disease, reacts, type="transformation", 
    correction_method="fdr", sorted=T, v_opt=NULL, full=FALSE) {
	# Create reference data
	cref = combn(1:ncol(normal), 2)
    normal = as.matrix(normal)
    disease = as.matrix(disease)
    
    if(type == "transformation") {
        normal = 1/normal
        disease = 1/disease
    }
    else if(type == "optimization") {
        M = cbind(normal,disease)
        S = get_stochiometry(reacts)
        if(is.null(v_opt)) v_opt = c(rep(0,ncol(S)-1),1)
        if (requireNamespace("foreach", quietly = TRUE)) {
            opt = foreach::"%dopar%"(foreach::foreach(i=1:ncol(M), .combine=cbind),  
                dba(v_opt, S, M[,i], lower=0, upper=1))
        }
        else {
            opt = lapply(1:ncol(M), function(i) 
                dba(v_opt, S, M[,i], lower=0, upper=1))
            opt = do.call(cbind, opt)
        }
        normal = opt[,1:ncol(normal)]
        disease = opt[,(ncol(normal)+1):ncol(M)]
    }
    else if(type != "raw") 
        stop("type must be either 'transformation' or 'optimization' :(") 
    
	res = single_hyp(normal[,1], normal[,1], reacts)
	res = res[,-ncol(res)]
	
	lfc_n = apply(cref, 2, function(idx) {
		h = single_hyp(normal[,idx[1]], normal[,idx[2]], reacts)
		return(h$log2_fold) })
	lfc_n = cbind(lfc_n, -lfc_n)
	
	# Create differential analysis
	ctreat = expand.grid(1:ncol(normal), 1:ncol(disease))
	
	lfc_d = apply(ctreat, 1, function(idx) {
		h = single_hyp(normal[,idx[1]], disease[,idx[2]], reacts)
		return(h$log2_fold) })
	
	# Generate statistics
	stats = lapply(1:nrow(res), function(i) {
		n_data = as.numeric(lfc_n[i,])
		d_data = as.numeric(lfc_d[i,])
		test = tryCatch({
				wilcox.test(x=d_data, y=n_data, conf.int=T)
			},
			error = function(e) {
				list(p.value=1, conf.int=rep(mean(n_data),2))
			},
			warning = function(w) {
				suppressWarnings(wilcox.test(x=d_data, y=n_data, conf.int=T))
			})
		
		return(data.frame(sd_normal=sd(n_data), sd_disease=sd(d_data),
				mean_log_fold=mean(d_data), ci_low=test$conf.int[1], 
				ci_high=test$conf.int[2], pval=test$p.value))
	})
	
	res = cbind(res, do.call(rbind, stats))
	reg = sapply(res$mean_log_fold, function(x) 
		if(x==0) "same" else if(x>0) "up" else "down")
	res$type = factor(reg)
	res$pval = p.adjust(res$pval, method=correction_method)
	if(sorted) res = res[order(res$pval, -abs(res$mean_log_fold)),]
	
	if(full) {
		lfc_n = data.frame(normal=lfc_n)
		lfc_d = data.frame(disease=lfc_d)
		res = list(hyp=res, lfc_normal=lfc_n, lfc_disease=lfc_d)
	}
	
	return(res)
}

