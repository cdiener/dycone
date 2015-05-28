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
	const_matrix = rcdd::d2q( -diag(ncol(s_matrix)) )
	const_b = rcdd::d2q(rep(0,ncol(s_matrix)))
	NC = rcdd::d2q( as.matrix( s_matrix%*%diag(mat) ) )
	b = rcdd::d2q(rep(0,nrow(s_matrix)))
	
	hp = rcdd::makeH(const_matrix, const_b, NC, b)
	hp = rcdd::addHeq(rep(1,ncol(s_matrix)), 1.0, hp)
	vrep = rcdd::scdd(hp)
	vp = vrep$output[,-(1:2)]
	
	basis = rcdd::q2d(t(vp))
	class(basis) = append("basis", class(basis))
	
	return( basis)
} 

get_stability = function(evs) {	
	evs[ abs(evs)<sqrt(.Machine$double.eps) ] = 0
	
	r = evs
	if ( all(Re(r)==0) && any(Im(r)!=0) ) return( "limit cycle" )
	if ( all(Re(r)==0) ) return ("steady state")
	if ( all(Re(r)<=0) ) return( "stable" )
	if ( all(Re(r)>=0) ) return ("unstable")
	else return ("undefined")
}

stability_analysis = function(basis, s_matrix, concs) {
	J = get_jacobian(s_matrix, concs)
	evs = NULL
	stab = NULL
	
	for( i in 1:ncol(basis) )
	{
		ev = eigen(s_matrix %*% diag(basis[,i]) %*% J)$values
		evs = rbind(evs, ev)
		stab = c(stab, get_stability(ev))
	}
	
	res = data.frame(what=stab, ev=evs, row.names=1:ncol(basis))
	names(res)[1:ncol(evs)+1] = paste0("ev",1:ncol(evs))
	
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

plot.basis = function(b, ...) {	
	pheatmap::pheatmap(b, border=NA, col=SINGLECOL(32), labels_row=1:nrow(b), 
		labels_col=1:ncol(b), ...)
}

plot.basis_diff = function(bd, ...) {
	pheatmap::pheatmap(bd, border=NA, col=DIVCOL(32), labels_row=1:nrow(bd), 
		labels_col=1:ncol(bd), ...)
}

plot_red = function(basis_list, col=c("seagreen", "yellow", "tomato")) {
	if( !("list" %in% class(basis_list)) ) stop("basis_list must be a list!")
	
	n = length(basis_list)
	all_basis = do.call(cbind, basis_list)
	pca = prcomp(t(all_basis))
	plot(NULL, xlim=c(-1,1), ylim=c(-1,1), xlab="PC 1", ylab="PC 2")
	pal = TRANSCOL(n)
	
	for( i in 1:n) {
		red = predict(pca, t(basis_list[[i]]))[,1:2]
		a = apply(red, 1, angle)
		red = red[order(a),]
		
		p = rbind(red,c(0,0))
		redundant_ids = rcdd::redundant( rcdd::makeV(rcdd::d2q(p)) )$redundant
		polygon(p[-redundant_ids,], border=NA, col=adjustcolor(pal[i], alpha.f=0.3))
		arrows(x0=0, y0=0, x1=red[,1], y1=red[,2], angle=15, length=0.05, col=pal[i])
	}
	energy = pca$sdev^2
	
	write(sprintf("Total standard deviation explained: %f%%.",
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

d = function(b1,b2) {
	if( all(dim(b1) == dim(b2)) ) res = b1-b2
	else if( nrow(b1) == nrow(b2) ) {
		bm = basis_map(b1,b2)
		res = b1[,bm$id1] - b2[,bm$id2]
	}
	else stop("Basis describe different reactions!")
	
	class(res) = append("basis_diff", class(res))
	
	return(res)
}

inside = function(x, s_matrix, v_terms, tol=sqrt(.Machine$double.eps)) {
	right = s_matrix%*%diag(v_terms)%*%x
	
	return( all(x>=0) & all(abs(right)<tol) )
}

eigenpathways = function(basis) {
	eps = svd(basis)$u
	eps[ abs(eps)<.Machine$double.eps ] = 0
	if( any(eps[,1]<0) ) eps = -eps
	
	return(eps)
}

hyp = function(b1, b2, reacts, log_fold_tol=0) {
	reacts = make_irreversible(reacts)
	if(!is.null(dim(b1))) b1 = rowMeans(b1)
	if(!is.null(dim(b2))) b2 = rowMeans(b2)
	
	zero_flux = b1==0 | b2==0
	b1 = b1[!zero_flux]
	b2 = b2[!zero_flux]
	logfold = log(b2,2) - log(b1,2)
	
	reg_id = which( abs(logfold)>=log_fold_tol )
	
	idx = id = r_s = kind = fac = NULL
	for( i in reg_id) {
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

all_diff = function(a,b) {
	combs = expand.grid(1:length(a), 1:length(b))
	diffs = apply(combs, 1, function(idx) a[idx[1]] - b[idx[2]])
	
	return(diffs)
}

multi_hyp = function(ref_list, treat_list, reacts, correction_method="fdr") {
	# Start by reducing the basis to row means
	ref = sapply(ref_list, rowMeans)
	treat = sapply(treat_list, rowMeans)
	
	# Create reference data
	cref = combn(1:ncol(ref), 2)
	print(cref)
	res = hyp(ref[,1], ref[,1], reacts)
	res = res[,-ncol(res)]
	
	lfc_ref = apply(cref, 2, function(idx) {
		h = hyp(ref[,idx[1]], ref[,idx[2]], reacts)
		return(h$log2_fold) })
	lfc_ref = cbind(lfc_ref, -lfc_ref)
	
	# Create differential analysis
	
	ctreat = expand.grid(1:ncol(ref), 1:ncol(treat))
	
	lfc_treat = apply(ctreat, 1, function(idx) {
		h = hyp(ref[,idx[1]], treat[,idx[2]], reacts)
		return(h$log2_fold) })
	
	# Generate statistics
	
	stats = lapply(1:nrow(res), function(i) {
		ref_data = as.numeric(lfc_ref[i,])
		treat_data = as.numeric(lfc_treat[i,])
		test = wilcox.test(x=treat_data, y=ref_data, conf.int=T)
		return(data.frame(mad_ref=mad(ref_data),
				mean_log_fold=mean(treat_data), 
				ci_low=test$conf.int[1], ci_high=test$conf.int[2],
				pval=test$p.value))
	})
	
	res = cbind(res, do.call(rbind, stats))
	reg = sapply(res$mean_log_fold, function(x) 
		if(x==0) "same" else if(x>0) "up" else "down")
	res$type = factor(reg)
	res$pval = p.adjust(res$pval, method=correction_method)
	res = res[order(res$pval),]
	
	return(res)
}
