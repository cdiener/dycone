#  ode.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

ion = function(x,idx) {
	return( if (idx %in% names(x)) x[idx] else NA )
}

get_sbml_params = function(sbml_file) {
	xml = XML::xmlParse(sbml_file)
	p_list = XML::xmlChildren(xml)[["sbml"]][["model"]][["listOfParameters"]]
	params = XML::xmlSApply(p_list, function(n) { 
		x = XML::xmlAttrs(n)
		return( c(ion(x,"id"), as.numeric(ion(x,"value"))) )
	})
	
	out = as.data.frame(t(params), row.names=FALSE)
	out$value = as.numeric(as.character(out$value))
	return( as.data.frame(t(params), row.names=FALSE) )
}

get_sbml_species = function(sbml_file) {
	xml = XML::xmlParse(sbml_file)
	s_list = XML::xmlChildren(xml)[["sbml"]][["model"]][["listOfSpecies"]]
	species = XML::xmlSApply(s_list, function(n) { 
		x = XML::xmlAttrs(n)
		return( c(ion(x,"id"), ion(x,"name"), ion(x, "initialConcentration")) )
	})
	
	out = as.data.frame(t(species), row.names=FALSE)
	out$initialConcentration = as.numeric(as.character(out$initialConcentration))
	return( out )
}

timecourse = function(x0, t, k, s_matrix) {
	f = function(t, y, p) { 
		list( s_matrix %*% diag(p) %*% get_ma_terms(s_matrix, y) )
	}
	
	sol = deSolve::lsoda(x0, t, f, k)
	
	return(sol)
}

plot.tc = function(tc) {
	melted = reshape2::melt(tc[,-1])
	melted = cbind(tc[,1], melted)
	names(melted) = c("t", "tid", "species", "conc")
	
	pl = ggplot2::ggplot(melted, ggplot2::aes(x=t, y=conc)) + ggplot2::geom_line() +
			facet_wrap(~species, scales="free_y") + ggplot2::theme_bw()
			
	return(pl)
}
