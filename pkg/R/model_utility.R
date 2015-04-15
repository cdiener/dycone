#  model_utility.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

# Helper functions in order to extract the model and construct basic data

str_trim = function(str) {
	return( gsub("^\\s+|\\s+$", "", str) )
}

#' Calculates the mass-action reaction rate
#'
#' @param substrates Stochiometry of the substrates
#' @param concs The concentrations used for the substrates
#' @return The mass action term \deqn{\prod_{i \in N^-} S_i^|N_i|}
mass_action = function(substrates, concs)
{
	sc = concs[substrates<0]
	nc = abs( substrates[substrates<0] )
	if(length(sc) == 0) return( 1 )
	
	return( prod(sc^nc) )
}

partial_deriv = function(i, substrates, concs)
{
	f = abs( min(substrates[i], 0) )
	substrates[i] = 0
	if( f==0 ) return (0)
	
	sc = concs[substrates<0]
	nc = abs( substrates[substrates<0] )
	
	if( length(sc) == 0 ) return( 0 )

	pd = prod(sc^nc)*f*concs[i]^(f-1)
	
	return( pd )
}

get_jacobian = function(s_matrix, concs)
{
	J = apply(s_matrix, 2, function(x) 
			sapply(seq_along(x), partial_deriv, substrates=x, concs=concs) )
	
	return( t(J) ) 
}

# gets the mass action terms \prod_j S_j^N(j,i)
get_ma_terms = function(s_matrix, concs)
{
	prods = apply( s_matrix, 2, mass_action, concs=concs)
	return( prods )
}

get_reaction_elems = function(reaction_str) {
	
	reversible = grepl("\\s*<.+\\s*", reaction_str)
	sides = strsplit(reaction_str, "\\s*<?\\s?=?-?\\s?>\\s*")[[1]]
	sides = strsplit(sides, "\\s*\\+\\s*")
	
	sub_pattern = "((\\d*)\\*|^)([^[:space:]]+)"
	subs = unlist( regmatches(sides[[1]], regexec(sub_pattern, sides[[1]])) )
	if (is.null(subs)) subs = c(NA, NA, "1", NA) 
	subs[subs==""] = "1"
	n_s = length(subs)
	
	if (length(sides)==1) prods = c(NA, NA, "1", NA)
	else {
		prods = unlist( regmatches(sides[[2]], regexec(sub_pattern, sides[[2]])) ) 
	}
	prods[prods==""] = "1"
	n_p = length(prods)
	
	return( list(S=subs[seq(4,n_s,4)], P=prods[seq(4,n_p,4)],
					N_S=as.numeric(subs[seq(3,n_s,4)]), 
					N_P=as.numeric(prods[seq(3,n_p,4)]), rev=reversible ) )
}

read_reactions = function(react_file) {
	reacts = read.csv(react_file)
	
	has_arrows = grepl("\\s*<?\\s?=?-?\\s?>\\s*", reacts[,1])
	if ( !all(has_arrows) ) {
		cat("Missing arrows in lines:")
		print(which(!has_arrows))
		stop("At least one of the lines has no reaction arrows!")
	}
	
	res = apply( reacts, 1, function(x) c( get_reaction_elems(x[1]), str_trim(x[-1]) ) )
	
	class(res) = append(class(res), "reactions")
	
	return( res )
}

build_reactions = function(s_matrix) {
	reacts = list()
	specs = rownames(s_matrix)
	for(i in 1:ncol(s_matrix)) {
		S = specs[ s_matrix[,i]<0 ]
		P = specs[ s_matrix[,i]>0 ]
		N_S = s_matrix[ s_matrix[,i]<0 , i ]
		N_P = s_matrix[ s_matrix[,i]>0 , i ]
		if (length(S) == 0) {
			S = NA
			N_S = 1
		}
		if (length(P) == 0) {
			P = NA
			N_P = 1
		}
		reacts[[i]] = list(S=S, P=P, N_S=abs(N_S), N_P=N_P, rev=FALSE, abbreviation=i)
	}
	
	class(reacts) = append(class(reacts), "reactions")
	
	return(reacts)
}

to_string = function(r, name=T) {
	id = paste0(r$abbreviation,": ")
	left = if (is.na(r$S[1])) "\u2205" 
			else paste(r$N_S, r$S, sep="*", collapse=" + ")
	right = if (is.na(r$P[1])) "\u2205" 
			else paste(r$N_P, r$P, sep="*", collapse=" + ")
	join = if (r$rev) "<=>" else "->"
	out = if(name) paste(id, left, join, right) else paste(left, join, right)
	return( out )
}

format.reactions = function(x) {
	r_str = sapply(x, to_string)
	
	return( paste(r_str, collapse="\n") ) 
}

print.reactions = function(x) {
	write( sprintf("Model has %d reactions (%d reversible)", 
			length(x), sum(sapply(x, function(a) a$rev))), file="" )
	write( format(x), file="" )
}

get_species = function(reacts) {
	if ( !("reactions" %in% class(reacts)) ) stop("Argument has wrong type!")
	
	species = unlist( lapply(reacts, function(x) c(x$S, x$P)) )
	species[is.na(species)] = "none"
	return( unique(species) )
}

get_stochiometry = function(reacts, reversible=FALSE, const="none") {
	if ( !("reactions" %in% class(reacts)) ) stop("Argument has wrong type!")
	
	species = get_species(reacts)
	n_r = if (reversible) length(reacts)
			else length(reacts) + sum(sapply(reacts, function(x) x$rev))
	
	N = matrix(0, nrow=length(species), ncol=n_r)
	rownames(N) = species
	i = 1
	for( r in reacts ) { 
		S = if (is.na(r$S)[1]) "none" else r$S
		P = if (is.na(r$P[1])) "none" else r$P
		N[S,i] = -r$N_S
		N[P,i] = r$N_P
		i = i+1
		if (!reversible && r$rev) {
			N[S,i] = r$N_S
			N[P,i] = -r$N_P
			i = i+1
		}
	}
	eliminate = rownames(N) %in% const
	N = N[!eliminate,]
	
	return(N)
}

as.reactions = function(s_matrix, reversible=F, r_names=NA) {
	if( class(rownames(s_matrix)) != "character" ) 
		stop("Not a valid stochiometric matrix (no rownames)!")
	
	if( length(r_names) == 1 && is.na(r_names) ) 
		r_names = paste0("r",1:ncol(s_matrix))
		
	species = rownames(s_matrix)
	reacts = vector("list", ncol(s_matrix))
	for( i in 1:ncol(s_matrix) ) {
		reacts[[i]]$S = species[s_matrix[,i]<0]
		if( length(reacts[[i]]$S) == 0 ) reacts[[i]]$S = NA
		reacts[[i]]$P = species[s_matrix[,i]>0]
		if( length(reacts[[i]]$P) == 0 ) reacts[[i]]$P = NA
		reacts[[i]]$N_S = -s_matrix[s_matrix[,i]<0,i]
		if( length(reacts[[i]]$N_S) == 0 ) reacts[[i]]$N_S = 1
		reacts[[i]]$N_P = s_matrix[s_matrix[,i]>0,i]
		if( length(reacts[[i]]$N_P) == 0 ) reacts[[i]]$N_P = 1
		
		if( length(reversible)==1 ) reacts[[i]]$rev = reversible
		else reacts[[i]]$rev = reversible[i]
		
		reacts[[i]]$name = reacts[[i]]$abbreviation = r_names[i]
	}
	
	class(reacts) = append("reactions", class(reacts))
	
	return(reacts)
}

make_irreversible = function(reacts) {
	r = list()
	for( i in 1:length(reacts) ) {
		rv = reacts[[i]]$rev
		reacts[[i]]$rev = FALSE
		
		r = append(r, list(reacts[[i]]))
		if ( rv ) {
			nr = reacts[[i]]
			nr$S = reacts[[i]]$P
			nr$P = reacts[[i]]$S
			nr$N_S = reacts[[i]]$N_P
			nr$N_P = reacts[[i]]$N_S
			r = append(r, list(nr)) 
		}
	}
	
	class(r) = append(class(r), "reactions")
	
	return(r)
}

#' Returns the number of different substrates a reaction has
#'
#' @param reacts A reactions object as returned by \code{\link{read_reactions}}
#' @return A numeric vector containing the number of substrates for each reaction
r_order = function(reacts) {
	return( sapply(reacts, function(x) length(x$S)) )
}

constant_flux = function(reacts) {
	return( sapply(reacts, function(x) any(is.na(x$S))) )
}

plot.reactions = function(x) {
	N = get_stochiometry(x, reversible=TRUE)
	if (requireNamespace("igraph", quietly = TRUE)) {
		g = igraph::graph.adjacency(N%*%t(N), weighted=TRUE, diag=FALSE)
		igraph::plot.igraph(g)
		
	} else {
		warning("igraph is not installed, Just showing the connectivity...")
		A = N%*%t(N)
		diag(A) = 0
		image(A!=0, col=c("white","black"))
	}
}

as.graph = function(reacts) {
	if ( !("reactions" %in% class(reacts)) ) stop("Argument has wrong type!")
	
	N = get_stochiometry(reacts, reversible=TRUE)
	adj = abs(N%*%t(N))
	
	if (requireNamespace("igraph", quietly = TRUE)) {
		return( igraph::graph.adjacency(adj, mode="max", 
				weighted=TRUE, diag=FALSE) )	
	} else {
		warning("igraph is not installed. Returning adjacency matrix...")
		return(adj)
	}
}

### Data documentation

#' Metabolic network of human red blood cell.
#'
#' Contains the reactions of the human red blood cell (erythrocite) metabolic 
#' model together with some randomly sample rates k.
#'
#' @format A list of reactions.
#' @source \url{http://journals.plos.org/plosone/article/asset?unique&id=info:doi/10.1371/journal.pone.0004967.s003}
"eryth"
