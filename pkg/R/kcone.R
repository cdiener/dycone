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
	all_basis = NULL
	for( b in basis_list ) all_basis = cbind(all_basis, b)
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
	write(sprintf("Total standard deviation explained: %f%%.",
				sum(pca$sdev[1:2])/sum(pca$sdev)*100), file="")
}

angle = function(x) {
	y = 0:1
	theta = acos( sum(x*y) / sqrt(sum(x * x)) ) 
	if(x[1]<0) theta = 2*pi - theta 
	return(theta)
}

basis_map = function(b1, b2) {
	closest = apply(b1, 2, function(x) {
							d = apply(b2, 2, function(y) dist(rbind(x,y)))
							i = which.min(d)
							return( c(i,d[i]) ) })
	
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

hyp = function(b1, b2, reacts, tol=5e-3) {
	reacts = make_irreversible(reacts)
	reg = rowMeans( d(b2,b1) ) 
	b1 = rowMeans(b1)
	b2 = rowMeans(b2)
	reg_id = which( abs(reg)>tol )
	
	idx = id = r_s = kind = fac = NULL
	for( i in reg_id) {
		r = reacts[[i]]
		idx = c(idx, i)
		id = c(id, r$abbreviation)
		r_s = c(r_s, to_string(r, name=F))
		kind = c(kind, if(reg[i]>0) "up" else "down")
		fac = c(fac,log(b2[i],2)-log(b1[i],2))
	}
	
	res = data.frame(idx=idx, name=id, reaction=r_s, type=kind, fold_change=fac)
	res = res[order(res$fold_change),]
	
	return(res)
}

