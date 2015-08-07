#  ode.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

#' Simulates the time course for a given model.
#'
#' @export
#' @keywords ODE, simulation
#' @param x0 Initial concentrations of metabolites.
#' @param t A vector of monotonously increasing times where the simulation
#'	is evaluated.
#' @param k The vector of kinetic constants used for simulation.
#' @param s_matrix The stochiometric matrix of the model.
#' @return A matrix containing the time in the first column and temporal
#'	concentrations of the metabolites in the following columns.
timecourse = function(x0, t, k, s_matrix) {
	f = function(t, y, p) { 
		list( s_matrix %*% diag(p) %*% get_ma_terms(s_matrix, y) )
	}
	
	sol = deSolve::lsoda(x0, t, f, k)
	class(sol) = append("dyconetc", class(sol))
	
	return(sol)
}
