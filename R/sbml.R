# sbml.R Copyright 2015 Christian Diener <ch.diener@gmail.com> MIT license.
# See LICENSE for more information.

# grep string to annotate the namespace
RE_SBML <- "sbml/level(\\d)/version(\\d)(/core|$)"
RE_FBC <- "sbml/level\\d/version\\d/fbc/version(\\d)"
RE_MATHML <- "Math/MathML"

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
        fbc=m_fbc[[1]][2]))
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
#'  xml_attr xml_text 
#' @export
sbml_params <- function(sbml_file) {
    attrs <- c("id", "name", "value")
    doc <- read_xml(sbml_file)
    ns <- clean_ns(xml_ns(doc))
    p_list <- xml_find_all(doc, 
        "./sbml:model/sbml:listOfParameters/sbml:parameter", 
        xml_ns(doc))
    params <- sapply(p_list, function(p) sapply(attrs, function(a) 
        xml_attr(p, a)))
    
    out <- as.data.frame(t(params), row.names = NULL)
    names(out) <- attrs
    out$value <- as.numeric(as.character(out$value))
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
    attrs <- c("id", "name", "initialamount", "initialconcentration")
    doc <- read_xml(sbml_file)
    ns <- clean_ns(xml_ns(doc))
    s_list <- xml_find_all(doc, 
        "./sbml:model/sbml:listOfSpecies/sbml:species", ns)
    species <- sapply(s_list, function(s) sapply(attrs, function(a) 
        xml_attr(s, a)))
    
    out <- as.data.frame(t(species), row.names = NULL)
    names(out) <- attrs
    num_cols <- grep("initial", names(out))
    
    for (i in num_cols) {
        out[, i] <- as.numeric(as.character(out[, i]))
    }
    return(out)
}

parse_reaction <- function(r_node, ns, v) {
    id <- r_node %>% xml_attr("id")
    name <- r_node %>% xml_attr("name")
    rev <- r_node %>% xml_attr("reversible")
    rev <- ifelse(rev == "true", TRUE, FALSE)
    
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
    
    # Get the notes
    notes <- r_node %>% xml_find_all("./sbml:notes", ns) %>% xml_text()
    
    # Get RDF annotations
    rdf_nodes <- r_node %>% xml_find_all("./sbml:annotation/rdf:RDF/rdf:Description/*", ns)
    rdf_qualifiers <- rdf_nodes %>% xml_name()
    rdf_links <- sapply(rdf_nodes, function(r) { 
        x <- r %>% xml_find_all("./rdf:Bag/rdf:li/@rdf:resource", ns) %>% xml_text()
        if (length(x) > 1) x <- sapply(x, paste0, sep=", ")
        x 
        })
    names(rdf_links) <- rdf_qualifiers
    
    r <- list(abbreviation=id, name=name, S=S, N_S=N_S, P=P, N_P=N_P, rev=rev)
    if (length(notes) > 0) r$notes <- notes
    if (length(rdf_links) > 0) r$rdf <- rdf_links
    return(r) 
}

#' Obtains the list of reactions from a SBML file.
#'
#' @param sbml_file A SBML file. This can be a file on the disk or the location
#'  to an online resource.
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
    doc <- read_xml(sbml_file)
    ns <- clean_ns(xml_ns(doc))
    v <- sbml_version(ns)
    
    r_nodes <- doc %>% xml_find_all("./sbml:model/sbml:listOfReactions/sbml:reaction", ns)
    
    if (progress) {
        cat("\n|")
        sapply(1:100, function(x) cat("-"))
        cat("|\n|=")
        trigger <- ceiling(length(r_nodes)/100)
        reacts <- lapply(1:length(r_nodes), function(i) {
            r <- parse_reaction(r_nodes[[i]], ns, v)
            if (i%%trigger == 0) cat("=")
            r
        })
        cat("|\n")
    }else reacts <- lapply(r_nodes, parse_reaction, ns=ns, v=v)
    class(reacts) <- append(class(reacts), "reactions")
    
    return(reacts)
}
