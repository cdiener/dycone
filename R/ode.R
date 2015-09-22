# ode.R Copyright 2015 Christian Diener <ch.diener@gmail.com> MIT license. See
# LICENSE for more information.

#' Simulates the time course for a given model.
#'
#' @export
#' @keywords ODE, simulation
#' @param x0 Initial concentrations of metabolites.
#' @param t A vector of monotonously increasing times where the simulation
#'  is evaluated.
#' @param k The vector of kinetic constants used for simulation.
#' @param S The stochiometric matrix of the model.
#' @return A matrix containing the time in the first column and temporal
#'  concentrations of the metabolites in the following columns.
#' @examples 
#' data(eryth)
#' S <- stochiometry(eryth)
#' x0 <- runif(nrow(S))
#' names(x0) <- rownames(S)
#' k <- runif(ncol(S))
#' t <- seq(0, 100, by=5)
#' timecourse(x0, t, k, S)
timecourse <- function(x0, t, k, S) {
    f <- function(t, y, p) {
        list(S %*% diag(p) %*% ma_terms(S, y))
    }
    
    sol <- deSolve::lsoda(x0, t, f, k)
    class(sol) <- append("dyconetc", class(sol))
    
    return(sol)
} 
