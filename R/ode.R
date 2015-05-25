#  ode.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

timecourse = function(x0, t, k, s_matrix) {
	f = function(t, y, p) { 
		list( s_matrix %*% diag(p) %*% get_ma_terms(s_matrix, y) )
	}
	
	sol = deSolve::lsoda(x0, t, f, k)
	class(sol) = append("dyconetc", class(sol))
	
	return(sol)
}

plot.dyconetc = function(tc) {
	melted = reshape2::melt(tc[,-1])
	melted = cbind(tc[,1], melted)
	names(melted) = c("t", "tid", "species", "conc")
	
	pl = ggplot2::ggplot(melted, ggplot2::aes(x=t, y=conc)) + ggplot2::geom_line() +
			facet_wrap(~species, scales="free_y") + ggplot2::theme_bw()
			
	return(pl)
}
