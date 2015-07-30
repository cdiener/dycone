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
