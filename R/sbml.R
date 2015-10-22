# sbml.R Copyright 2015 Christian Diener <ch.diener@gmail.com> MIT license.
# See LICENSE for more information.

# grep string to annotate the namespace
RE_SBML <- "sbml/level(\\d)/version(\\d)(/core|$)"
RE_FBC <- "sbml/level\\d/version\\d/fbc/version(\\d)"
RE_MATHML <- "Math/MathML"

get_xml <- function(obj) {
    if ("xml_node" %in% class(obj)) return(obj)
    else return(read_xml(obj))
}

# Get the sbml version
sbml_version <- function(ns) {
    sbml_ns <- ns["sbml"]
    fbc_ns <- ns["fbc"]
    re <- regexec(RE_SBML, sbml_ns)
    m_sbml <- regmatches(sbml_ns, re)
    if ("fbc" %in% names(ns)) {
        re <- regexec(RE_FBC, fbc_ns)
        m_fbc <- regmatches(fbc_ns, re)
    } else m_fbc <- list(c(NA, NA))
    return(c(level=as.numeric(m_sbml[[1]][2]), version=as.numeric(m_sbml[[1]][3]),
        fbc=as.numeric(m_fbc[[1]][2])))
}

# Clean up the namespace
clean_ns <- function(ns) {
    sbml_idx <- grep(RE_SBML, ns)
    fbc_idx <- grep(RE_FBC, ns)
    mml_idx <- grep(RE_MATHML, ns)
    
    if (length(sbml_idx) != 1) stop("No or several SBML namespaces found.")
    else names(ns)[sbml_idx] <- "sbml"
    if (length(fbc_idx) == 1) names(ns)[fbc_idx] <- "fbc"
    if (length(mml_idx) == 1) names(ns)[mml_idx] <- "mathml"
    
    return(ns)
}

#' Obtains the list of parameters from a SBML model
#'
#' @param sbml_file A SBML file. This can be a file on the disk or the location
#'  to an online resource.
#' @return A data frame containing the obtained data. Some entries might be NA.
#' @examples
#' # requires internet connection
#' m_url <- "http://www.ebi.ac.uk/compneur-srv/biomodels-main/download?mid=MODEL1103210001"
#' sbml_params(m_url)
#'
#' @importFrom xml2 read_xml read_html xml_find_one xml_find_all xml_ns 
#'  xml_attr xml_attrs xml_name xml_text 
#' @export
sbml_params <- function(sbml_file) {
    attrs <- c("id", "name", "value")
    doc <- get_xml(sbml_file)
    ns <- clean_ns(xml_ns(doc))
    p_list <- doc %>% xml_find_all("./sbml:model/sbml:listOfParameters/sbml:parameter", ns)
    params <- sapply(p_list, function(p) sapply(attrs, function(a) 
        xml_attr(p, a)))
    
    out <- as.data.frame(t(params), row.names = NULL)
    if (ncol(out) == 3) {
        names(out) <- attrs
        out$value <- as.numeric(as.character(out$value))
        if (any(duplicated(out$id))) stop("There are duplicated parameters IDs!")
    } else return(NULL)
    return(out)
}

#' Obtains the list of species from a SBML model
#'
#' @param sbml_file A SBML file. This can be a file on the disk or the location
#'  to an online resource.
#' @return A data frame containing the obtained data. Some entries might be NA.
#' @examples
#' # requires internet connection
#' m_url <- "http://www.ebi.ac.uk/compneur-srv/biomodels-main/download?mid=MODEL1103210001"
#' sbml_species(m_url)
#'
#' @export
sbml_species <- function(sbml_file) {
    doc <- get_xml(sbml_file)
    ns <- clean_ns(xml_ns(doc))
    v <- sbml_version(ns)
    
    if(is.na(v[3])) {
        attrs <- c("id", "name", "initialAmount", "initialConcentration", 
            "boundaryCondition", "charge")
    } else {
        attrs <- c("id", "name", "initialAmount", "initialConcentration", 
            "boundaryCondition", "fbc:charge", "fbc:chemicalFormula")
    }
    
    s_list <- doc %>% xml_find_all("./sbml:model/sbml:listOfSpecies/sbml:species", ns)
    species <- sapply(s_list, function(s) sapply(attrs, function(a) 
        xml_attr(s, a, ns)))
    
    out <- as.data.frame(t(species), row.names = NULL)
    names(out) <- attrs
    num_cols <- c(3, 4, 6)
    
    for (i in num_cols) {
        out[, i] <- as.numeric(as.character(out[, i]))
    }
    out$boundaryCondition <- ifelse(out$boundaryCondition == "true", TRUE, FALSE)
    out$boundaryCondition[is.na(out$boundaryCondition)] <- FALSE
    return(out)
}

parse_reaction <- function(r_node, ns, v) {
    id <- r_node %>% xml_attr("id")
    name <- r_node %>% xml_attr("name")
    rev <- r_node %>% xml_attr("reversible")
    rev <- ifelse(rev == "true", TRUE, FALSE)
    if (is.na(rev) && v[1] < 3) rev <- TRUE
    
    # Get the substrates and products
    subs <- r_node %>% xml_find_all("./sbml:listOfReactants/
        sbml:speciesReference", ns) %>% xml_attrs()
    S <- sapply(subs, "[", "species")
    N_S <- sapply(subs, function(x) as.numeric(x["stoichiometry"]))
    if (length(S) == 0) {
        S <-  NA
        N_S <- 1
    }
    if (any(is.na(N_S))) {
        if(v[1] > 2) {
            stop(sprintf("At least one substrate in reaction %s has no stoichiometry!", id))
        } else N_S[is.na(N_S)] <- 1
    }
    names(S) <- NULL
    names(N_S) <- NULL
    
    prods <- subs <- r_node %>% xml_find_all("./sbml:listOfProducts/
        sbml:speciesReference", ns) %>% xml_attrs()
    P <- sapply(prods, "[", "species")
    N_P <- sapply(prods, function(x) as.numeric(x["stoichiometry"]))
    if (length(P) == 0) {
        P <-  NA
        N_P <- 1
    }
    if (any(is.na(N_P))) {
        if(v[1] > 2) {
            stop(sprintf("At least one product in reaction %s has no stoichiometry!", id))
        } else N_P[is.na(N_P)] <- 1
    }
    names(P) <- NULL
    names(N_P) <- NULL
    
    # Get bounds
    if (!is.na(v[3])) {
        lower <- r_node %>% xml_attr("fbc:lowerFluxBound", ns)
        upper <- r_node %>% xml_attr("fbc:upperFluxBound", ns)
    } else {
        lower <- r_node %>% xml_find_all(
        "./sbml:kineticLaw/sbml:listOfParameters/sbml:parameter[@id='LOWER_BOUND']", ns) %>%
        xml_attr("value")
        upper <- r_node %>% xml_find_all(
        "./sbml:kineticLaw/sbml:listOfParameters/sbml:parameter[@id='UPPER_BOUND']", ns) %>%
        xml_attr("value")
    }
    if (length(lower) == 0) lower <- NA
    if (length(upper) == 0) upper <- NA
    
    # Get the notes
    notes <- r_node %>% xml_find_all("./sbml:notes", ns) %>% xml_text()
    if (length(notes) == 0) notes <- NA
    
    # Get RDF annotations
    rdf_nodes <- r_node %>% xml_find_all("./sbml:annotation/rdf:RDF/rdf:Description/*", ns)
    rdf_qualifiers <- rdf_nodes %>% xml_name()
    rdf_links <- sapply(rdf_nodes, function(r) { 
        x <- r %>% xml_find_all("./rdf:Bag/rdf:li/@rdf:resource", ns) %>% xml_text()
        if (length(x) > 1) x <- sapply(x, paste0, sep=", ")
        x 
        })
    if (length(rdf_links) == 0) rdf_links <- NA
    else names(rdf_links) <- rdf_qualifiers
    
    r <- list(abbreviation=id, name=name, S=S, N_S=N_S, P=P, N_P=N_P, rev=rev, 
        notes=notes, rdf=rdf_links, lower=lower, upper=upper)
    return(r) 
}

#' Obtains the list of reactions from a SBML file.
#'
#' @param sbml_file A SBML file. This can be a file on the disk or the location
#'  to an online resource.
#' @param progress Whether a progress indication should be shown.
#' @return A list with additional class "reactions" describing the reactions in
#'  the model. RDF annotations, notes and flux bounds are saved along with the 
#'  reactions
#' @examples
#' # requires internet connection
#' m_url <- "http://www.ebi.ac.uk/compneur-srv/biomodels-main/download?mid=MODEL1103210001"
#' sbml_reactions(m_url)
#'
#' @export
sbml_reactions <- function(sbml_file, progress=TRUE) {
    doc <- get_xml(sbml_file)
    ns <- clean_ns(xml_ns(doc))
    v <- sbml_version(ns)
    
    r_nodes <- doc %>% xml_find_all("./sbml:model/sbml:listOfReactions/sbml:reaction", ns)
    
    if (progress) {
        n <- length(r_nodes)
        trigger <- ceiling(length(r_nodes)/100)
        reacts <- lapply(1:n, function(i) {
            cat("                                                           \r")
            cat(sprintf("Parsing reaction %d/%d", i, n))
            parse_reaction(r_nodes[[i]], ns, v)
        })
        cat(" -- Done.\n")
    } else reacts <- lapply(r_nodes, parse_reaction, ns=ns, v=v)
    
    class(reacts) <- append(class(reacts), "reactions")
    return(reacts)
}

#' Main driver to parse an SBML model.
#'
#' This function wraps around \code{sbml_*} functions but also performs some 
#' basic sanity checks and adjustments which are:
#' \itemize{
#' \item{Gives information about the used level and version and whether the FBC
#'  package is used.}
#' \item{Checks whether there are species that do not appear in reactions (warning).}
#' \item{Checks whether all species references in reactions are defined in the
#'  species list (error).}
#' \item{Maps lower and upper flux bounds to its respective global parameters.}
#' }
#'
#' @seealso \code{\link{sbml_species}}, \code{\link{sbml_params}} and 
#'  \code{\link{sbml_reactions}} to read only parts of an SBML model.
#' @param sbml_file A SBML file. This can be a file on the disk or the location
#'  to an online resource.
#' @return A list containing three elements:
#'  \describe{
#'  \item{species}{The species (metabolites or molecules) in the model in a data
#'  frame. Boundary metabolites are removed automatically.}
#'  \item{params}{The global parameters defined for the model.}
#'  \item{reactions}{The reactions in the model as an "reactions" object.} 
#'  }
#' @examples
#' # requires internet connection
#' m_url <- "http://www.ebi.ac.uk/compneur-srv/biomodels-main/download?mid=MODEL1103210001"
#' read_sbml(m_url)
#'
#' @export
read_sbml <- function(sbml_file) {
    doc <- get_xml(sbml_file)
    ns <- clean_ns(xml_ns(doc))
    v <- sbml_version(ns)
    write(sprintf("Reading SBML level %d version %d %s FBC.", v[1], v[2],
        ifelse(is.na(v[3]), "without", "with")), file="")
    
    specs <- sbml_species(doc)
    bc <- specs$boundaryCondition
    if (sum(bc) > 0) warning("Boundary species will be ignored!")
    specs <- specs[!specs$boundaryCondition, ]
    
    reacts <- sbml_reactions(doc, progress=TRUE)
    
    # Check for species sanity
    r_specs <- species(reacts)
    if (all(is.na(specs$id))) spec_list <- specs$name
    else spec_list <- specs$id 
    not_in_rs <- which(!(spec_list %in% r_specs))
    not_in_slist <- which(!(r_specs %in% spec_list))
    if (length(not_in_rs) > 0) warning(
        paste("The following species are defined but not used in reactions:",
        paste(spec_list[not_in_rs], collapse=", ")))
    if (length(not_in_rs) > 0) stop(
        paste("The following species used in reactions are undefined:",
        paste(r_specs[not_in_slist], collapse=", ")))
    
    # Map flux bounds
    p <- sbml_params(doc)
    if(!is.null(p)) {
        cat("Mapping flux bounds to global parameters")
        fixed <- sapply(1:length(reacts), function(i) {
            r <- reacts[[i]]
            fixed <- FALSE
            if(is.character(r$lower)) {
                if (!(r$lower %in% p$id)) stop(paste0("There is no global parameter ",
                    r$lower, "."))
                reacts[[i]]$lower <<- p$value[p$id == r$lower]
                fixed <- TRUE
            }
            if(is.character(r$upper)) {
                if (!(r$lower %in% p$id)) stop(paste0("There is no global parameter ",
                    r$lower, "."))
                reacts[[i]]$upper <<- p$value[p$id == r$upper] 
                fixed <- TRUE
            }
            fixed
        })
        write(paste(" --", sum(fixed), "reactions mapped."), file="")
    } else write("No global parameters defined.", file="")
    
    return(list(species=specs, reactions=reacts, params=p))
}
