# model_utility.R Copyright 2015 Christian Diener <ch.diener@gmail.com> MIT
# license. See LICENSE for more information.

# Helper functions in order to extract the model and construct basic data

### General helper functions
str_conv = function(str) {
    val = tryCatch(as.numeric(str), warning = function(w) {
        if (is.null(grep(",", str))) 
            return(gsub("^\\s+|\\s+$", "", str))
        v = gsub("^\\s+|\\s+$", "", strsplit(str, ",[[:space:]]*")[[1]])
        v = tryCatch(as.numeric(v), warning = function(w) v)
        return(v)
    })
    
    return(val)
}

order_by = function(x, y) {
    o = sapply(x, function(x) which(x == y)[1])
    return(o)
}

#' %c% is a caching operator. It is used after an
#' expression, detects all assigned variables and saves
#' them into a file. If the file given exists the variables
#' are read from the cache file and not calculated again.
#' In order to delete the cache simply delete the associated
#' file.
#'
#' @export
#' @keywords operator, caching
#' @param ex The R expression to be cached. Use brackets for 
#'\tmultiline statements.\t
#' @param filename The filename of the cache file. 
#' @examples
#' # Is executed only if the file does not exist, otherwise
#' # samples a and b are read from 'samples.Rd'
#' { a <- rnorm(1e5) } %c% 'samples.Rd'
"%c%" = function(ex, filename) {
    if (!file.exists(filename)) {
        e = new.env()
        eval(match.call()$ex, envir = e)
        save(list = ls(e), file = filename, envir = e)
    }
    
    load(filename, envir = parent.frame(4))
}

#' Calculates a mass-action reaction rate
#'
#' \code{mass_action} is mostly used inetrnally to get mass action rates for a 
#' single column of the stochiometric matrix. It does not check whether the 
#' ordering of substrates in the stochiometry and concs is corrects, so please
#' make sure that the indices in \code{substrates} and \code{concs} coincide.  
#'
#' @seealso \code{\link{stochiometry}} to generate the stochiometric
#'\tmatrix from a reaction list. \code{\link{ma_terms}} for a wrapper that
#'\tcalculates all mass-action terms for a stochiometric matrix and also checks 
#'\tfor correct ordering of the substrates.
#'\t\code{\link{ma_terms}} to calculate all mass-action terms at once.
#' @export
#' @param substrates Stochiometry of the substrates. Must have the same length
#'\tas concs.
#' @param concs The concentrations used for the substrates
#' @return The mass action term \deqn{\prod_{i \in N^-} S_i^|N_i|}
#' @examples
#' data(eryth)
#' S = stochiometry(eryth)
#' ma1 = mass_action(S[,1], runif(nrow(S)))
mass_action = function(substrates, concs) {
    if (length(substrates) != length(concs)) 
        stop("substrates and concs must have the same length!")
    
    sc = concs[substrates < 0]
    nc = abs(substrates[substrates < 0])
    if (length(sc) == 0) 
        return(1)
    
    return(prod(sc^nc))
}

#' Calculates derivatives of mass-action kinetics.
#'
#' \code{deriv_ma} does not check whether the 
#' ordering of substrates in the stochiometry and concs is corrects, so please
#' make sure that the indices in \code{substrates} and \code{concs} coincide.
#'
#' @seealso \code{\link{jacobian}} to calculate the Jacobian matrix.
#' @export
#' @param i The index of the substrate used as the differentiation variable.
#' @param substrates The stochiometry of the subbstrates (column in the
#'\tstochiometric matrix. Must be the same length as concs.
#' @param concs Concentrations of the substrates.
#' @examples
#' data(eryth)
#' S = stochiometry(eryth)
#' dma1_dS3 = deriv_ma(3, S[,1], runif(nrow(S)))
deriv_ma = function(i, substrates, concs) {
    if (length(substrates) != length(concs)) 
        stop("substrates and concs must have the same length!")
    
    f = abs(min(substrates[i], 0))
    substrates[i] = 0
    if (f == 0) 
        return(0)
    
    sc = concs[substrates < 0]
    nc = abs(substrates[substrates < 0])
    
    # if( length(sc) == 0 ) return( 0 )
    
    pd = prod(sc^nc) * f * concs[i]^(f - 1)
    
    return(pd)
}

#' Calculates the Jacobian matrix.
#'
#' @seealso \code{\link{stochiometry}} to calculate the stochiometric matrix.
#' @export
#' @keywords jacobian, stochiometry
#' @param s_matrix The stochiometric matrix of the model.
#' @param concs The concentrations of the substrates. Must be a named vector with
#'\tnames being the metabolite names as used in \code{rownames(s_matrix)}.
#' @param deriv_func The function calculating the derivative. Must be of the 
#'\tform func(idx, subs, concs) where idx is the index of the differentiation 
#'\tvariable, subs the stochiometry of the substrates and concs the concentration
#'\tof substrates. See \code{\link{mass_action}} for an example.
#' @return The Jacobian matrix of dimension n_s x n_r where entry (i,j) denotes
#'\tthe derivative dv[j]/dS[i].
#' @examples
#' data(eryth)
#' S = stochiometry(eryth)
#' concs = runif(nrow(S))
#' names(concs) = rownames(S)
#' J = jacobian(S, concs)
jacobian = function(s_matrix, concs, deriv_func = deriv_ma) {
    if (is.null(names(concs))) 
        stop("Concentration vector must have names!") else {
        concs = concs[rownames(s_matrix)]
        prods = apply(s_matrix, 2, mass_action, concs = concs)
    }
    
    J = apply(s_matrix, 2, function(x) sapply(seq_along(x), deriv_func, substrates = x, 
        concs = concs))
    
    return(t(J))
}

#' Calculates all mass action terms \deqn{\prod_{i \in N^-} S_i^|N_i|}
#'
#' @seealso \code{\link{stochiometry}} to calculate the stochiometric matrix.
#' @export
#' @keywords kinetics, mass-action
#' @param s_matrix The stochiometric matrix.
#' @param concs The concentrations of the substrates. Can be one of two
#'\t\itemize{
#'\t\t\item A named vector containing the concentratiions for all substrates, 
#'\t\twhere its names come from the same set as the rownames of S.
#'\t\t\item A data frame with one column names 'name' containing the 
#'\t\tnames and all other columns containing multiple concentration entries.
#'\t}
#' @return Either vector of length n_r containing the mass-action terms, if
#'\tconcs was a vector, or a matrix with n_r rows containing the mass-action
#'\tin its columns.
#' @examples
#' data(eryth)
#' S = stochiometry(eryth)
#' concs = runif(nrow(S))
#' names(concs) = rownames(S)
#' mats = ma_terms(S, concs)
ma_terms = function(s_matrix, concs) {
    prods = NULL
    
    if (is.vector(concs)) {
        if (is.null(names(concs))) 
            stop("Concentration vector must have names!") else {
            concs = concs[rownames(s_matrix)]
            prods = apply(s_matrix, 2, mass_action, concs = concs)
        }
    } else if (is.data.frame(concs)) {
        if ("name" %in% names(concs)) {
            o = order_by(rownames(s_matrix), concs$name)
            concs = concs[o, ]
            name_idx = which(names(concs) == "name")
            prods = apply(concs[, -name_idx], 2, function(co) apply(s_matrix, 2, 
                mass_action, concs = co))
        } else stop("Concentration data frame must have a column named 'name'!")
    }
    
    return(prods)
}

# Internal function to get the stochiometry from a string representation of the
# reaction
get_reaction_elems = function(reaction_str) {
    
    reversible = grepl("\\s*<.+\\s*", reaction_str)
    sides = strsplit(reaction_str, "\\s*<?\\s?=?-?\\s?>\\s*")[[1]]
    sides = strsplit(sides, "\\s*\\+\\s*")
    
    sub_pattern = "((\\d*\\.*\\d*)\\*|^\\s*)([^[:space:]]+)"
    subs = unlist(regmatches(sides[[1]], regexec(sub_pattern, sides[[1]])))
    if (is.null(subs)) 
        subs = c(NA, NA, "1", NA)
    subs[subs == ""] = "1"
    n_s = length(subs)
    
    if (length(sides) == 1) 
        prods = c(NA, NA, "1", NA) else {
        prods = unlist(regmatches(sides[[2]], regexec(sub_pattern, sides[[2]])))
    }
    prods[prods == ""] = "1"
    n_p = length(prods)
    
    return(list(S = subs[seq(4, n_s, 4)], P = prods[seq(4, n_p, 4)], N_S = as.numeric(subs[seq(3, 
        n_s, 4)]), N_P = as.numeric(prods[seq(3, n_p, 4)]), rev = reversible))
}

#' Reads a list of reactions with optional annotations from a file.
#' 
#' @export
#' @param react_file A csv file containing the reactions, where the first
#'\tcolumn denotes the reaction string and additional columns contain
#'  annotations. A single annotation can contain several values as a 
#'\tcomma separated string. The csv file should be saved with all entries
#'\tquoted.
#' @examples
#' r_str = 'reaction,abbreviation,numbers\nA -> B,blub,'1,2,3'\nB <=>, bla, 3'
#' r = read_reactions(textConnection(r_str))
read_reactions = function(react_file) {
    reacts = read.csv(react_file, stringsAsFactors = FALSE)
    
    has_arrows = grepl("\\s*<?\\s?=?-?\\s?>\\s*", reacts[, 1])
    if (!all(has_arrows)) {
        stop(sprintf("The following reactions are missing reacion arrows: %s", paste(which(!has_arrows), 
            collapse = ", ")))
    }
    
    res = apply(reacts, 1, function(x) c(get_reaction_elems(x[1]), lapply(x[-1], 
        str_conv)))
    
    class(res) = append(class(res), "reactions")
    
    return(res)
}

# Helper function to convert a reaction entry to a clean string version
to_string = function(r, name = T) {
    id = paste0(r$abbreviation, ": ")
    left = if (is.na(r$S[1])) 
        "∅" else paste(r$N_S, r$S, sep = "*", collapse = " + ")
    right = if (is.na(r$P[1])) 
        "∅" else paste(r$N_P, r$P, sep = "*", collapse = " + ")
    join = if (r$rev) 
        "<=>" else "->"
    out = if (name) 
        paste(id, left, join, right) else paste(left, join, right)
    return(out)
}

# Formats reactions to a nice string representation
format.reactions = function(x) {
    r_str = sapply(x, to_string)
    
    return(paste(r_str, collapse = "\n"))
}

#' Prints a nice description of the reactions.
#'
#' @seealso \code{\link{as.reactions}} to convert a stochiometric matrix
#'\tto a reaction list.
#' @export
#' @param x The reaction list to be output.
#' @param ... Additional arguments passed on to write.
#' @examples
#' data(eryth)
#' print(eryth)
print.reactions = function(x, ...) {
    write(sprintf("Model has %d reactions (%d reversible)", length(x), sum(sapply(x, 
        function(a) a$rev))), file = "")
    write(format(x), file = "", ...)
}

#' Gets all species/metabolites from a reaction list.
#'
#' @export
#' @keywords susbtrates, stochiometry
#' @param reacts A reaction list.
#' @return A vector containing all unique species in the model.
#' @examples
#' data(eryth)
#' print(species(eryth))
species = function(reacts) {
    if (!("reactions" %in% class(reacts))) 
        stop("Argument has wrong type!")
    
    species = unlist(lapply(reacts, function(x) c(x$S, x$P)))
    species[is.na(species)] = "none"
    return(unique(species))
}

#' Calculates the stochiometric matrix for a list of reactions.
#' 
#' @seealso \code{\link{read_reactions}} to read a reaction list from a file.
#' @export
#' @keywords stochiometry
#' @param reacts The reaction list to be used.
#' @param reversible Whether the stochiometric matrix can include reversible 
#'\treactions.
#' @param const A vector of species names that are assumed to be invariant and
#'\twill be dropped from the stochiometric matrix. The default is not to drop
#'\tany species.
#' @return The stochiometric matrix with dimension n_s x n_r.
#' @examples
#' data(eryth)
#' S = stochiometry(eryth)
stochiometry = function(reacts, reversible = FALSE, const = NULL) {
    if (!("reactions" %in% class(reacts))) 
        stop("Argument has wrong type!")
    
    spec = species(reacts)
    n_r = if (reversible) 
        length(reacts) else length(reacts) + sum(sapply(reacts, function(x) x$rev))
    
    N = matrix(0, nrow = length(spec), ncol = n_r)
    rownames(N) = spec
    i = 1
    for (r in reacts) {
        S = if (is.na(r$S)[1]) 
            "none" else r$S
        P = if (is.na(r$P[1])) 
            "none" else r$P
        N[S, i] = -r$N_S
        N[P, i] = r$N_P
        i = i + 1
        if (!reversible && r$rev) {
            N[S, i] = r$N_S
            N[P, i] = -r$N_P
            i = i + 1
        }
    }
    eliminate = rownames(N) %in% c("none", const)
    N = N[!eliminate, ]
    
    return(N)
}

#' Converts a stochiometric matrix to a reaction list.
#'
#' @export
#' @keywords stochiometry, reactions
#' @param s_matrix The stochiometric matrix to be converted.
#' @param reversible Marks reversible reactions. FALSE denotes that
#'\tall reactions are irreversible. Otherwise a boolean vector of length
#'\tn_r defining the reversibility for each reaction.
#' @param r_names If NA reaction names are generated as r1,...,rn. Can be a a vector
#'  of length n_r denoting names for the reactions.
#' @return A reaction list containing the reactions from S.
#' @examples
#' S = matrix(c(-1,1,1,-1), nrow=2)
#' rownames(S) <- c('A', 'B')
#' print(as.reactions(S))
as.reactions = function(s_matrix, reversible = F, r_names = NA) {
    if (class(rownames(s_matrix)) != "character") 
        stop("Not a valid stochiometric matrix (no rownames)!")
    
    if (length(r_names) == 1 && is.na(r_names)) 
        r_names = paste0("r", 1:ncol(s_matrix))
    
    species = rownames(s_matrix)
    reacts = vector("list", ncol(s_matrix))
    for (i in 1:ncol(s_matrix)) {
        reacts[[i]]$S = species[s_matrix[, i] < 0]
        if (length(reacts[[i]]$S) == 0) 
            reacts[[i]]$S = NA
        reacts[[i]]$P = species[s_matrix[, i] > 0]
        if (length(reacts[[i]]$P) == 0) 
            reacts[[i]]$P = NA
        reacts[[i]]$N_S = -s_matrix[s_matrix[, i] < 0, i]
        if (length(reacts[[i]]$N_S) == 0) 
            reacts[[i]]$N_S = 1
        reacts[[i]]$N_P = s_matrix[s_matrix[, i] > 0, i]
        if (length(reacts[[i]]$N_P) == 0) 
            reacts[[i]]$N_P = 1
        
        if (length(reversible) == 1) 
            reacts[[i]]$rev = reversible else reacts[[i]]$rev = reversible[i]
        
        reacts[[i]]$name = reacts[[i]]$abbreviation = r_names[i]
    }
    
    class(reacts) = append("reactions", class(reacts))
    
    return(reacts)
}

#' Converts a given reaction list into an irreversible one, splitting up
#' all reversible ones into two irreversible reactions.
#' 
#' @export
#' @keywords reactions
#' @param reacts The reaction list.
#' @return A reaction list containing all reactions in reacts with reversible
#'  reactions being split in two.
#' @examples
#' data(eryth)
#' print(make_irreversible(eryth))
make_irreversible = function(reacts) {
    r = list()
    for (i in 1:length(reacts)) {
        rv = reacts[[i]]$rev
        reacts[[i]]$rev = FALSE
        
        r = append(r, list(reacts[[i]]))
        if (rv) {
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
#' @export
#' @param r The reaction list.
#' @param field Name of the property to be obtained.
#' @return A data frame mapping the property to reaction indices.
#' @examples
#' data(eryth)
#' print(rp(eryth, 'name'))
rp = function(r, field = "KEGG_enzyme") {
    if (length(field) != 1) 
        stop("field must be exactly one string!")
    prop = lapply(1:length(r), function(i) {
        ri = r[[i]]
        if (all(is.na(ri[[field]]))) 
            return(NULL)
        data.frame(r_idx = i, x = ri[[field]], stringsAsFactors = FALSE)
    })
    prop = do.call(rbind, prop)
    names(prop)[2] = field
    
    return(prop)
}

#' Returns the number of different substrates a reaction has
#'
#' @export
#' @param reacts A reactions object as returned by \code{\link{read_reactions}}
#' @return A numeric vector containing the number of reactants for each reaction
#' @examples
#' data(eryth)
#' print(r_order(eryth))
r_order = function(reacts) {
    return(sapply(reacts, function(x) sum(x$N_S) - is.na(x$S[1])))
}

#' Identifies reactions with a constant flux, thus reactions with no
#' substrates.
#'
#' @seealso \code{\link{r_order}} to get the reaction order.
#' @export
#' @keywords reactions
#' @param reacts The reaction list.
#' @return A boolean vector identifying reactions with constant flux.
#' @examples
#' data(eryth)
#' which(constant_flux(eryth))
constant_flux = function(reacts) {
    return(sapply(reacts, function(x) any(is.na(x$S))))
}

#' Plots the reactions as a graph.
#' 
#' @seealso \code{\link{print.reactions}} for a nicely formatted representation.
#' @export
#' @keywords plot, reactions
#' @param x A reaction list.
#' @param ... Additional arguments passed to plot.igraph.
#' @examples
#' data(eryth)
#' plot(eryth)
plot.reactions = function(x, ...) {
    N = stochiometry(x, reversible = TRUE)
    
    if (requireNamespace("igraph", quietly = TRUE)) {
        g = igraph::graph.adjacency(N %*% t(N), weighted = TRUE, diag = FALSE)
        igraph::plot.igraph(g, layout = igraph::layout.circle, vertex.size = 10, 
            edge.arrow.size = 0.5, ...)
        
    } else {
        warning("igraph is not installed, Just showing the connectivity...")
        A = N %*% t(N)
        diag(A) = 0
        image(A != 0, col = c("white", "black"))
    }
}

#' Converts a reaction list into a graph object.
#'
#' @export
#' @keywords reactions, graph
#' @param reacts A reaction list.
#' @param ... other arguments passed to specific methods
#' @return An igraph object representing with the species being the nodes and
#'  reactions denoting edges.
#' @examples
#' data(eryth)
#' print(as.graph(eryth))
as.graph = function(reacts) {
    if (!("reactions" %in% class(reacts))) 
        stop("Argument has wrong type!")
    
    N = stochiometry(reacts, reversible = TRUE)
    adj = abs(N %*% t(N))
    
    if (requireNamespace("igraph", quietly = TRUE)) {
        return(igraph::graph.adjacency(adj, mode = "max", weighted = TRUE, diag = FALSE))
    } else {
        warning("igraph is not installed. Returning adjacency matrix...")
        return(adj)
    }
} 
