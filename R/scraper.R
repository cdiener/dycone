#  scraper.py
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

HMDB_SEARCH = "http://www.hmdb.ca/unearth/q?query=%s&searcher=metabolites"
HMDB_XML = "http://www.hmdb.ca/metabolites/%s.xml"
KEGG_REST = "http://rest.kegg.jp/link/hsa/ec:%s"

#' Obtains the list of parameters from a SBML model
#'
#' @param sbml_file A SBML file. This can be a file on the disk or the location
#' 	to an online resource.
#' @return A data frame containing the obtained data. Some entries might be NA.
#' @export
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
#' @export
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
	val = gsub(",", "", as.character(val))
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
#' @export
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
    
    empty_or_text = function(x) {if(is.null(x)) "" else rvest::xml_text(x)}    
    
	vals = sapply(nodes, function(n) sapply(tags, function(ta) {
		rvest::xml_node(n, ta) %>% empty_or_text()
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
#' @param hmids The HMDB ID or a vector of IDs for the metabolites.
#' @param add A data frame with as many rows as entries in hmid containing additional
#'	information that will be added to the result (for instance names or ids).
#' @return The scraped data set as a data frame or NULL if no concentrations were
#'	found.
#' @export
hmdb_concentration = function(hmids, add=NULL) {
	out = NULL
	cat("\n")
	for(i in 1:length(hmids)) {
		id = hmids[i]
		cat(sprintf("\rScraping %d/%d...",i,length(hmids)))
		hm_xml = rvest::xml(sprintf(HMDB_XML,id))
		
		hm_entries = hm_xml %>% 
			rvest::xml_nodes("normal_concentrations concentration") %>% hmdb_parse()
		
		kegg_id = hm_xml %>% rvest::xml_node("kegg_id") %>% rvest::xml_text()
		name = hm_xml %>% rvest::xml_node("name") %>% rvest::xml_text()
		
		if(!is.null(hm_entries)) {	
			hm_entries = cbind(kegg_id, id, name, hm_entries)
			names(hm_entries)[c(1:3,ncol(hm_entries))] = c("keggid", "hmdbid","name", "pmid")
			if(is.data.frame(add)) hm_entries = cbind(hm_entries, add[i,])
			out = rbind(out, hm_entries)
		}
	}
	
	return(out)
}

# Primitive function for group patching missing data
group_patch = function(x) {
	miss = t(apply(x, 1, is.na))
	
	n_m = rowSums(miss)
	fixable = which(n_m > 0 & n_m < ncol(x))
	for(i in fixable) {
		nas = is.na(x[i,])
		x[i,nas] = mean(as.numeric(x[i,!nas]))
	}
	
	return(x)
}

# Primitive function for reference patching missing data
ref_patch = function(x, ref) {
	miss = t(apply(x[,-1], 1, is.na))
	
	need_fix = which(rowSums(miss) > 0)
	for(i in need_fix) {
		id = as.character(x[i,1])
		nas = is.na(x[i,])
		nas[1] = FALSE
		ref_data = as.numeric(unlist(ref[ref[,1] == id, -1]))
		x[i, nas] = mean(ref_data, na.rm=T)
	} 
	
	return(x)
}

#' Patches holes in a data set by using group means and reference data.
#'
#' Biological data is only rarely complete, however k-cone analysis requires
#' fully populated steady state concentration data. This functions helps to patch
#' holes in your data. It assumes that the measurements data set is a data frame
#' where the id column uses a shared identifier as the reference data set and
#' the data columns specifify control and treatment measurements.
#' The optional reference data set is further used to fill up holes which cannot
#' be patched by group means. 
#'
#' @param measurements A data frame with one id column, and at least one normal 
#' 	and one treatment column
#' @param id An index or column name specifying the ids
#' @param normal A vector of indices or column names specifying the control 
#'	group
#' @param treatment A vector of indices or column names specifying the treatment 
#'	group 
#' @param ref_data An optional data frame of reference data 2 or 3 columns, 
#'	where the first column denotes the ids, the second the normal reference data 
#'	and the (optional) third column the treatment reference data.
#' @return The original measurement data frame with a maximum number of missing 
#' data being patched.
#' @export
patch = function(measurements, id, normal, treatment, ref_data=NULL) {
	out = measurements
	data_cols = c(normal,treatment)
	
	out[,normal] = group_patch(out[,normal])		# patch grouped by normal
	out[,treatment] = group_patch(out[,treatment])	# patch grouped by treatment
	
	out[,data_cols] = group_patch(out[,data_cols])	# patch between groups
	
	# patch normal group by reference
	if(!is.null(ref_data)) {
		out[,data_cols] = ref_patch(out[,c(id,data_cols)], ref_data)[,-1]
	}
	
	return(out)
}

