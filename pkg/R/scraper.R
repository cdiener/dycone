#  scraper.py
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

HMDB_SEARCH = "http://www.hmdb.ca/unearth/q?query=%s&searcher=metabolites"
HMDB_XML = "http://www.hmdb.ca/metabolites/%s.xml"

#' Obtains the list of parameters from a SBML model
#'
#' @param sbml_file A SBML file. This can be a file on the disk or the location
#' 	to an online resource.
#' @return A data frame containing the obtained data. Some entries might be NA.
sbml_params = function(sbml_file) {
	attrs = c("id", "name", "value")
	doc = rvest::xml(sbml_file)
	p_list = doc %>% rvest::xml_nodes("sbml model listofparameters parameter")
	params = sapply(p_list, function(p) 
		sapply(attrs, function(a) rvest::xml_attr(p,a)) )
	
	out = as.data.frame(t(params), row.names=NULL)
	names(out) = attrs
	out$value = as.numeric(as.character(out$value))
	return( out )
}

#' Obtains the list of species from a SBML model
#'
#' @param sbml_file A SBML file. This can be a file on the disk or the location
#' 	to an online resource.
#' @return A data frame containing the obtained data. Some entries might be NA.
sbml_species = function(sbml_file) {
	attrs = c("id", "name", "initialamount", "initialconcentration")
	doc = rvest::xml(sbml_file)
	s_list = doc %>% rvest::xml_nodes("sbml model listofspecies species")
	species = sapply(s_list, function(s) 
		sapply(attrs, function(a) rvest::xml_attr(s,a)) )
	
	out = as.data.frame(t(species), row.names=NULL)
	names(out) = attrs
	num_cols = grep("initial", names(out))
	
	for( i in num_cols ) {
		out[,i] = as.numeric(as.character(out[,i]))
	}
	return( out )
}

#' Parses concentration values from HMDB.
#'
#' @param val The value as it appears in HMDB.
#' @return The numeric concentration value.
parse_conc = function(val) {
	val = as.character(val)
	val = strsplit(val, " ")[[1]][1]
	if(length(grep("-",val))==1) {
		vals = as.numeric(strsplit(val,"-")[[1]])
		val = mean(vals)
	}
	return(as.numeric(val))
}

#' Finds IDs of compounds from HMDB corresponding to a given search term.
#'
#' @param search_term Can be any text search term or a different ID (for instance
#' 	a KEGG ID
#' @return A single HMDB ID or list of putative IDs.
find_hmdb = function(search_term) {
	hmids = rvest::html(sprintf(HMDB_SEARCH, RCurl::curlEscape(search_term))) %>%
			rvest::html_nodes(".result-link .btn-card") %>% rvest::html_text()
			
	return(hmids)
}

#' Parses basic data from a set of HMDB concentration nodes.
#'
#' @param nodes A list of XML nodes.
#' @return The parsed data set as a data frame.
hmdb_parse = function(nodes) {
	tags = c("biofluid", "concentration_value", "concentration_units", 
		"subject_age", "subject_sex", "subject_condition", 
		"references reference pubmed_id")
	vals = sapply(nodes, function(n) sapply(tags, function(ta) {
		rvest::xml_node(n, ta) %>% {if(is.null(.)) "" else rvest::xml_text(.)}
		}) )
	
	vals = t(tolower(vals))
	out = as.data.frame(vals, row.names=NULL)
	
	if(ncol(out) == 0) return(NULL)
	
	names(out) = tags
	out$concentration_value = sapply(out$concentration_value, parse_conc)
	
	return(out)
}

#' Scrapes measured metabolite concentration values for a given HMDB ID.
#'
#' @param hmid The HMDB ID of the metabolite.
#' @return The scraped data set as a data frame or NULL if no concentrations were
#'	found.
hmdb_concentration = function(hmid) {
	hm_xml = rvest::xml(sprintf(HMDB_XML,hmid))
	
	hm_entries = hm_xml %>% 
		rvest::xml_nodes("normal_concentrations concentration") %>% hmdb_parse()
	
	kegg_id = hm_xml %>% rvest::xml_node("kegg_id") %>% rvest::xml_text()
	name = hm_xml %>% rvest::xml_node("name") %>% rvest::xml_text()
	
	hm_entries = cbind(kegg_id, hmid, name, hm_entries)
	
	if(is.null(hm_entries)) return(NULL)
	else {	
		hm_entries = cbind(kegg_id, hmid, name, hm_entries)
		names(hm_entries)[c(1:3,ncol(hm_entries))] = c("keggid", "hmdbid","name", "pmid")
		return(hm_entries)
	}
}
