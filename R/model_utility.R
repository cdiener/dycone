#  model_utility.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

# Helper functions in order to extract the model and construct basic data

### General helper functions
str_conv = function(str) {
	val = tryCatch(as.numeric(str), warning = function(w) {
		if( is.null(grep(",",str)) )
			return(gsub("^\\s+|\\s+$", "", str))
		v = gsub("^\\s+|\\s+$", "", strsplit(str,",[[:space:]]*")[[1]])
        v = tryCatch(as.numeric(v), warning = function(w) v)
        return(v)
	})

	return(val)		
}

order_by = function(x, y) {
	o = sapply(x, function(x) which(x==y)[1])
	return(o)
}

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
'%c%' = function(ex, filename) {
	if(!file.exists(filename)) {
		e = new.env()
		eval(match.call()$ex, envir=e)
		save(list=ls(e), file=filename, envir=e)
	}
	
	load(filename, env=parent.frame(4)) 
}

#' Calculates the mass-action reaction rate
#'
#' @export
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

deriv_ma = function(i, substrates, concs)
{
	f = abs( min(substrates[i], 0) )
	substrates[i] = 0
	if( f==0 ) return (0)
	
	sc = concs[substrates<0]
	nc = abs( substrates[substrates<0] )
	
	# if( length(sc) == 0 ) return( 0 )

	pd = prod(sc^nc)*f*concs[i]^(f-1)
	
	return( pd )
}

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
get_jacobian = function(s_matrix, concs, deriv_func = deriv_ma)
{
	J = apply(s_matrix, 2, function(x) 
			sapply(seq_along(x), deriv_func, substrates=x, concs=concs) )
	
	return( t(J) ) 
}

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
get_ma_terms = function(s_matrix, concs)
{
	prods = NULL
	
	if(is.vector(concs)) {
		if(is.null(names(concs))) stop("Concentration vector must have names!")
		else {
			concs = concs[rownames(s_matrix)]
			prods = apply( s_matrix, 2, mass_action, concs=concs)
		}
	}
	else if(is.data.frame(concs)) {
		if("name" %in% names(concs)) {
			o = order_by(concs$name, rownames(S))
			concs = concs[o,]
			name_idx = which(names(concs) == "name")
			prods = apply(concs[,-name_idx], 2, function(co) 
						apply( s_matrix, 2, mass_action, concs=co))
		}
		else stop("Concentration data frame must have a column named 'name'!")
	}
	
	return( prods )
}

get_reaction_elems = function(reaction_str) {
	
	reversible = grepl("\\s*<.+\\s*", reaction_str)
	sides = strsplit(reaction_str, "\\s*<?\\s?=?-?\\s?>\\s*")[[1]]
	sides = strsplit(sides, "\\s*\\+\\s*")
	
	sub_pattern = "((\\d*\\.*\\d*)\\*|^\\s*)([^[:space:]]+)"
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
	reacts = read.csv(react_file, stringsAsFactors=FALSE)
	
	has_arrows = grepl("\\s*<?\\s?=?-?\\s?>\\s*", reacts[,1])
	if ( !all(has_arrows) ) {
		stop(sprintf("The following reactions are missing reacion arrows: %s",
            paste(which(!has_arrows), collapse=", ")))
	}
	
	res = apply(reacts, 1, function(x) 
		c( get_reaction_elems(x[1]), lapply(x[-1], str_conv) ))
	
	class(res) = append(class(res), "reactions")
	
	return( res )
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

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
get_species = function(reacts) {
	if ( !("reactions" %in% class(reacts)) ) stop("Argument has wrong type!")
	
	species = unlist( lapply(reacts, function(x) c(x$S, x$P)) )
	species[is.na(species)] = "none"
	return( unique(species) )
}

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
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

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
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

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
make_irreversible = function(reacts) {
	r = list()
	for( i in 1:length(reacts) ) {
		rv = reacts[[i]]$rev
		reacts[[i]]$rev = FALSE
		
		r = append(r, list(reacts[[i]]))
		if ( rv ) {
			br = reacts[[i]]
			br$S = reacts[[i]]$P
			br$P = reacts[[i]]$S
			br$N_S = reacts[[i]]$N_P
			br$N_P = reacts[[i]]$N_S
			r = append(r, list(br)) 
		}
	}
	
	class(r) = append(class(r), "reactions")
	
	return(r)
}

#' Obtains properties from a reaction list
#' 
#' @param r The reaction list.
#' @param field Name of the property to be obtained.
#' @return A data.frame mapping the property to reaction indices.
#' @export
rp = function(r, field="KEGG_enzyme") {
	prop = lapply(1:length(r), function(i) {
        ri = r[[i]]
        if(all(is.na(ri[[field]]))) return(NULL)
        data.frame(r_idx=i, x=ri[[field]], stringsAsFactors=F)
    })
    prop = do.call(rbind,prop)
    names(prop)[2] = field
    
	return(prop)
}

#' Returns the number of different substrates a reaction has
#'
#' @param reacts A reactions object as returned by \code{\link{read_reactions}}
#' @return A numeric vector containing the number of substrates for each reaction
#' @export
r_order = function(reacts) {
	return( sapply(reacts, function(x) length(x$S) - is.na(x$S[1])) )
}

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
constant_flux = function(reacts) {
	return( sapply(reacts, function(x) any(is.na(x$S))) )
}

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
plot.reactions = function(x) {
	N = get_stochiometry(x, reversible=TRUE)

	if (requireNamespace("igraph", quietly = TRUE)) {
		g = igraph::graph.adjacency(N%*%t(N), weighted=TRUE, diag=FALSE)
		igraph::plot.igraph(g, layout=igraph::layout.circle, vertex.size=10, 
            edge.arrow.size=0.5)
		
	} else {
		warning("igraph is not installed, Just showing the connectivity...")
		A = N%*%t(N)
		diag(A) = 0
		image(A!=0, col=c("white","black"))
	}
}

#' TODO: change me >:(
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(df, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(df)}
#'    \item \code{ggplot()}
#'   }
#' 
#' @seealso \url{http://had.co.nz/ggplot2}
#'  \code{\link{geom_segment}} for a more general approach
#' @export
#' @keywords internal
#' @param data default data set
#' @param ... other arguments passed to specific methods
#' @return The stuff :O
#' @examples
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
#'                  y = rnorm(30))
#' # Compute sample mean and standard deviation in each group
#' library(plyr)
#' ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
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
