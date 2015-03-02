#  kcone.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

# Functions to analyze the k-cone

get_kcone_basis = function(s_matrix, v_terms) {	
	mat = v_terms
	SM = s_matrix %*% diag(mat)
	
	return( MASS::Null( t(SM) ) )
}

get_polytope_basis = function(s_matrix, v_terms) {
	mat = v_terms
	const_matrix = rcdd::d2q( -diag(ncol(s_matrix)) )
	const_b = rcdd::d2q(rep(0,ncol(s_matrix)))
	NC = rcdd::d2q( as.matrix( s_matrix%*%diag(mat) ) )
	b = rcdd::d2q(rep(0,nrow(s_matrix)))
	
	hp = rcdd::makeH(const_matrix, const_b, NC, b)
	hp = rcdd::addHeq(rep(1,ncol(s_matrix)), 1.0, hp)
	vp = rcdd::scdd(hp)$output[,-(1:2)]
	
	return( rcdd::q2d(t(vp)) )
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

plot_basis = function(b_matrix) {	
	#rownames(b_matrix) = reactions
	bdf = reshape2::melt( b_matrix )
	names(bdf) = c("reaction", "idx", "value")
	
	pl = ggplot2::ggplot(bdf, ggplot2::aes(x=idx, y=reaction, fill=value)) + 
		ggplot2::geom_raster() + ggplot2::theme_bw() +
		ggplot2::scale_x_continuous(expand=c(0,0)) + 
		ggplot2::scale_y_continuous(expand=c(0,0)) +
		if (min(b_matrix)>=0) ggplot2::scale_fill_continuous(low="white", high="darkred")
		else ggplot2::scale_fill_gradient2()

	return(pl)
}
