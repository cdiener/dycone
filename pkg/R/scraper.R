#  scraper.py
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

HMDB_SEARCH = "http://www.hmdb.ca/unearth/q?query=%s&searcher=metabolites"
HMDB_XML = "http://www.hmdb.ca/metabolites/%s.xml"

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

parse_conc = function(val) {
	val = as.character(val)
	val = strsplit(val, " ")[[1]][1]
	if(length(grep("-",val))==1) {
		vals = as.numeric(strsplit(val,"-")[[1]])
		val = mean(vals)
	}
	return(as.numeric(val))
}

hmdb_parse = function(nodes) {
	tags = c("biofluid", "concentration_value", "concentration_units", 
		"subject_age", "subject_sex", "subject_condition", 
		"references reference pubmed_id")
	vals = sapply(nodes, function(n) sapply(tags, function(ta) {
		rvest::xml_node(n, ta) %>% {if(is.null(.)) "" else rvest::xml_text(.)}
		}) )
	
	vals = t(tolower(vals))
	out = as.data.frame(vals, row.names=NULL)
	names(out) = tags
	out$concentration_value = sapply(out$concentration_value, parse_conc)
	
	return(out)
}

hmdb_concentration = function(search_term) {
	hmid = rvest::html(sprintf(HMDB_SEARCH, search_term)) %>%
			rvest::html_node(".result-link .btn-card") %>% rvest::html_text()
	
	hm_xml = rvest::xml(sprintf(HMDB_XML,hmid))
	
	hm_entries = hm_xml %>% 
		rvest::xml_nodes("normal_concentrations concentration") %>% hmdb_parse()
	
	kegg_id = hm_xml %>% rvest::xml_node("kegg_id") %>% rvest::xml_text()
	name = hm_xml %>% rvest::xml_node("name") %>% rvest::xml_text()
				
	hm_entries = cbind(kegg_id, hmid, name, hm_entries)
	names(hm_entries)[c(1:3,ncol(hm_entries))] = c("keggid", "hmdbid","name", "pmid")
	
	return(hm_entries)
}
