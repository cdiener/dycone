#  kcone.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

# Functions to analyze the k-cone

SINGLECOL = colorRampPalette(c("white","darkred"))
DIVCOL = colorRampPalette(c("blue", "white", "tomato"))
TRANSCOL = colorRampPalette(c("seagreen", "darkgoldenrod1", "tomato"))

get_kcone_basis = function(s_matrix, v_terms) {	
	mat = v_terms
	SM = s_matrix %*% diag(mat)
	
	basis = MASS::Null( t(SM) )
	class(basis) = append("basis", class(basis))
	
	return( basis )
}

get_polytope_basis = function(s_matrix, v_terms) {
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
	
	basis = rcdd::q2d(t(vp))
	class(basis) = append("basis", class(basis))
	
	return( basis)
} 

save_ine = function(s_matrix, v_terms) {
	
}

get_stability = function(evs) {	
	evs[ abs(evs)<sqrt(.Machine$double.eps) ] = 0
	
	r = evs
	if ( all(Re(r)==0) && any(Im(r)!=0) ) return( "limit cycle" )
	if ( all(Re(r)==0) ) return ("zero jacobian")
	if ( all(Re(r)<=0) ) return( "stable" )
	else return ("unstable")
}

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

get_fluxes = function(basis, S, reactions, concs) {
	mat = get_ma_terms(S, concs)
	imp = rowMeans(basis)*mat
	names(imp) = 1:length(imp)

	rn = NULL
	for(r in reactions) {
		rn = rbind(rn, c(r$abbreviation, r$name))
		if( r$rev ) rn = rbind(rn, c(r$abbreviation, r$name))
	}

	imp = data.frame( flux=imp, id=rn[,1], name=rn[,2] )
	imp = imp[order(imp$flux, decreasing=T),]
	
	return(imp)
}

occupation = function(basis) {
	if( !("basis" %in% class(basis)) ) stop("object is not a basis!")
	
	return(apply(basis, 1, function(r) 
		sum(r>sqrt(.Machine$double.eps))/length(r)))
}

plot.basis = function(b, ...) {	
	pheatmap::pheatmap(b, border=NA, col=SINGLECOL(32), labels_row=1:nrow(b), 
		labels_col=1:ncol(b), ...)
}

plot_red = function(basis_list, arrows=TRUE, col=NULL, n_cl=NULL) {
	if( !("list" %in% class(basis_list)) ) stop("basis_list must be a list!")
	
	b_rows = sapply(basis_list, nrow)
    b_cols = sapply(basis_list, ncol)
    
    if(sd(b_rows)>0) stop("All basis need to have the same dimension!")
    if(nrow(basis_list[[1]])<2) stop("Basis must be at least 2-dimensional!")
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
    
    if(max(b_cols)>1e3) {
        write("Basis are very large. Reducing by clustering...", file="")
        if(is.null(n_cl)) n_cl = sqrt(mean(b_cols))
        cl = lapply(red, function(b) 
            kmeans(b,centers=n_cl, iter.max=100, nstart=1))
        
        dists = sapply(cl, function(x) mean(sqrt(x$withinss/x$size)))
        write(sprintf("Mean in-cluster distance: %g.",mean(dists)), file="")
        
        red = lapply(cl, function(x) x$centers)
    }
    
	rs = apply(do.call(rbind,red), 2, range)
	plot(NULL, xlim=rs[,1], ylim=rs[,2], xlab="PC 1", ylab="PC 2")
	
	if(is.null(col)) pal = TRANSCOL(n)
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

angle = function(x) {
	y = 0:1
	theta = acos( sum(x*y) / sqrt(sum(x * x)) ) 
	if(x[1]<0) theta = 2*pi - theta 
	return(theta)
}

basis_map = function(b1, b2) {
	n1 = ncol(b1)
	n2 = ncol(b2)
	progress = F
	cur = 0
	
	if(n1*n2>1e5) {
		write("Basis are pretty large. This might take a long time!", file="")
		progress = T
	}
	closest = apply(b1, 2, function(x) {
							d = apply(b2, 2, function(y) dist(rbind(x,y)))
							i = which.min(d)
							if(progress) {
								cur <<- cur+1
								cat("\r") 
								cat(sprintf("Finished %.2f%%...",cur/n1*100))
							}
							return( c(i,d[i]) )
							
							 })
	cat("\n")
	
	return( data.frame(id1=1:ncol(b1), id2=closest[1,], d=closest[2,]) )
}

inside = function(x, s_matrix, v_terms, tol=sqrt(.Machine$double.eps)) {
	right = s_matrix%*%diag(v_terms)%*%x
	
	if(is.null(dim(x))) return(all(right>=0) & all(abs(right)<tol))
	return( apply(right, 2, function(r) all(r>=0) & all(abs(r)<tol)) )
}

eigendynamics = function(basis, n=1) {
	if(is.null(dim(basis))) return(basis)
	
	eps = svd(basis)$u
	eps[ abs(eps)<.Machine$double.eps ] = 0
	if( any(eps[,1]<0) ) eps = -eps
	
	eps = apply(eps,2,function(x) x/sum(x))
	
	return(eps[,1:n])
}

r_means = function(x) {
	if(!is.null(dim(x))) x = rowMeans(x)
	
	return(x)
}

hyp = function(b1, b2, reacts, log_fold_tol=0) {
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
#' @param ref_list A list of basis for the reference/control group.
#' @param treat_list A list of basis for the treatment/disease group.
#' @param reacts A reaction list.
#' @param correction_method A correction method for the multiple test p-values.
#'	Takes the same arguments as the method argument in p.adjust.
#' @param full If TRUE also returns the individual log fold changes along with
#'	the differential regulation data.
#' @return If full is TRUE returns a list of generated hypothesis and the individual
#' 	log fold changes between all reference basis and between reference and treatments.
#'	If full is FALSE only returns the generated hypothesis. 
multi_hyp = function(normal, disease, reacts, correction_method="fdr", full=FALSE) {
	# Create reference data
	cref = combn(1:ncol(normal), 2)
	res = hyp(normal[,1], normal[,1], reacts)
	res = res[,-ncol(res)]
	
	lfc_n = apply(cref, 2, function(idx) {
		h = hyp(normal[,idx[1]], normal[,idx[2]], reacts)
		return(h$log2_fold) })
	lfc_n = cbind(lfc_n, -lfc_n)
	
	# Create differential analysis
	ctreat = expand.grid(1:ncol(normal), 1:ncol(disease))
	
	lfc_d = apply(ctreat, 1, function(idx) {
		h = hyp(normal[,idx[1]], disease[,idx[2]], reacts)
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
	res = res[order(res$pval, -abs(res$mean_log_fold)),]
	
	if(full) {
		lfc_n = data.frame(normal=lfc_n)
		lfc_d = data.frame(disease=lfc_d)
		res = list(hyp=res, lfc_normal=lfc_n, lfc_disease=lfc_d)
	}
	
	return(res)
}
