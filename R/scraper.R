# scraper.R Copyright 2015 Christian Diener <ch.diener@gmail.com> MIT license.
# See LICENSE for more information.

HMDB_SEARCH <- "http://www.hmdb.ca/unearth/q?query=%s&searcher=metabolites"
HMDB_XML <- "http://www.hmdb.ca/metabolites/%s.xml"
KEGG_REST <- "http://rest.kegg.jp/link/hsa/ec:%s"


#' Extracts mean concentrations from HMDB results using a selection of biofluids
#' in a distinct order.
#'
#' This is useful if one prefers concentrations from one biofluid, but will also
#' except other biofluids if no reported measurements can be found for the
#' primary one. \code{priority_mean} will return the mean of the first biofluid
#' for which measuremnts can be found.
#'
#' @export
#' @seealso \code{\link{hmdb_concentration}} to parse the concentrations used
#'  in this function.
#' @param d A data frame as returned by hmdb_concentration
#' @param biofluids Vector of valid biofluids to be checked in the order
#'  they appear.
#' @return The mean
#' @examples
#' # requires internet connection
#' concs <- hmdb_concentration('HMDB00124')
#' priority_mean(concs) # has measurements for cytoplasm so returns those
priority_mean <- function(d, biofluids = c("cellular cytoplasm", "blood")) {
    if (is.null(d))
        return(NA)

    m <- NULL

    for (bf in biofluids) {
        m <- d$concentration_value[d$biofluid == bf]
        if (length(m) > 0 && !all(is.na(m)))
            break
    }

    return(mean(m, na.rm = TRUE))
}

#' Wrapper for grep that returns NA if nothing was found. Also splits up
#' comma sperated arguments.
#'
#' @export
#' @param ids ID or comma separated IDs to be searched.
#' @param x Character in which to search.
#' @return The result of grep or NA if not found.
#' @examples
#' grep_id('A, B', 'B, D')
grep_id <- function(ids, x) {
    res <- unlist(sapply(str_conv(ids), grep, x = x))

    if (length(res) > 0) return(res) else return(NA)
}

#' Parses concentration values from HMDB.
#'
#' @param val The value as it appears in HMDB.
#' @return The numeric concentration value.
parse_conc <- function(val) {
    val <- gsub(",", "", as.character(val))
    val <- strsplit(val, " ")[[1]][1]
    if (length(grep("-", val)) == 1) {
        vals <- as.numeric(strsplit(val, "-")[[1]])
        val <- mean(vals)
    }
    return(as.numeric(val))
}

#' Finds IDs of compounds from HMDB corresponding to a given search term.
#'
#' @param search_term Can be any text search term or a different ID (for instance
#'  a KEGG ID
#' @return A single HMDB ID or list of putative IDs.
#' @examples
#' # requires internet connection
#' find_hmdb("pyruvate")
#'
#' @importFrom RCurl curlEscape
#' @export
find_hmdb <- function(search_term) {
    hmids <- read_html(sprintf(HMDB_SEARCH, curlEscape(search_term))) %>%
        xml_find_all("//div[@class='result-link']/a[@class='btn-card']") %>%
        xml_text()

    return(hmids)
}

#' Parses basic data from a set of HMDB concentration nodes.
#'
#' @param nodes A list of XML nodes.
#' @return The parsed data set as a data frame.
hmdb_parse_concentration <- function(nodes) {
    tags <- c("biofluid", "concentration_value", "concentration_units",
              "subject_age", "subject_sex", "subject_condition",
              "references/reference/pubmed_id")

    vals <- sapply(nodes, function(n) sapply(tags, function(ta) {
        tryCatch(n %>% xml_find_first(paste0("./", ta)) %>% xml_text(),
            error = function(e) return(""))
    }))

    vals <- t(tolower(vals))
    out <- as.data.frame(vals, row.names = NULL)

    if (ncol(out) == 0)
        return(NULL)

    names(out) <- tags
    out$concentration_value <- sapply(out$concentration_value, parse_conc)

    return(out)
}

#' Parses the XML entry for a single HMDB id
#'
#' @param id The id for the entry
#' @return The entry as a data frame
hmdb_parse_xml <- function(id) {
    hm_xml <- read_xml(sprintf(HMDB_XML, id))

    entry <- hm_xml %>%
        xml_find_all("./normal_concentrations/concentration") %>%
        hmdb_parse_concentration()

    kegg_id <- hm_xml %>% xml_find_first("./kegg_id") %>%
        xml_text()
    name <- hm_xml %>% xml_find_first("./name") %>% xml_text()
    if (!is.null(entry)) {
        entry <- cbind(kegg_id, id, name, entry, row.names = NULL)
        names(entry)[c(1:3, ncol(entry))] <- c("keggid",
            "hmdbid", "hmdb_name", "pmid")
    }

    return(entry)
}

#' Scrapes measured metabolite concentration values for a given HMDB ID.
#'
#' @param hmids The HMDB ID or a vector of IDs for the metabolites.
#' @param add A data frame with as many rows as entries in hmid containing additional
#'  information that will be added to the result (for instance names or ids).
#' @return The scraped data set as a data frame or NULL if no concentrations were
#'  found.
#' @examples
#' # requires internet connection
#' hmdb_concentration('HMDB00124')
#'
#' @export
#' @importFrom xml2 xml_find_first xml_find_all
hmdb_concentration <- function(hmids, add = NULL) {
    out <- NULL
    cat("\n")
    for (i in 1:length(hmids)) {
        ids <- str_conv(hmids[i])
        hm_entries <- NULL
        kegg_id <- NULL
        name <- NULL
        for (id in ids) {
            cat(sprintf("\rScraping %d/%d...", i, length(hmids)))
            entry <- hmdb_parse_xml(id)
            if (is.data.frame(add) && !is.null(entry))
                entry <- cbind(entry, add[i, ], row.names = NULL)
            hm_entries <- rbind(hm_entries, entry)
        }
        out <- rbind(out, hm_entries)
    }
    cat("\n")

    return(out)
}

# Primitive function for group patching missing data
group_patch <- function(x) {
    x <- as.matrix(x)
    miss <- t(apply(x, 1, is.na))

    n_m <- rowSums(miss)
    fixable <- which(n_m > 0 & n_m < ncol(x))
    for (i in fixable) {
        nas <- is.na(x[i, ])
        x[i, nas] <- mean(as.numeric(x[i, !nas]))
    }

    return(x)
}

# Primitive function for reference patching missing data
ref_patch <- function(x, ref) {
    y <- as.matrix(x[, -1])
    miss <- t(apply(y, 1, is.na))

    need_fix <- which(rowSums(miss) > 0)
    for (i in need_fix) {
        id <- as.character(x[i, 1])
        nas <- is.na(y[i, ])
        ref_data <- as.numeric(unlist(ref[ref[, 1] == id, -1]))
        y[i, nas] <- mean(ref_data, na.rm = TRUE)
    }
    x[, -1] <- y

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
#'  and one treatment column
#' @param id An index or column name specifying the ids
#' @param normal A vector of indices or column names specifying the control
#'  group
#' @param treatment A vector of indices or column names specifying the treatment
#'  group
#' @param ref_data An optional data frame of reference data 2 or 3 columns,
#'  where the first column denotes the ids, the second the normal reference data
#'  and the (optional) third column the treatment reference data.
#' @return The original measurement data frame with a maximum number of missing
#'  data being patched.
#' @export
#' @examples
#' x <- data.frame(id = 1:3, normal = c(NA, 1, 2), disease = c(2, NA, 3))
#' patch(x, 1, 2, 3)
patch <- function(measurements, id, normal, treatment, ref_data = NULL) {
    out <- measurements
    data_cols <- c(normal, treatment)

    # patch grouped by normal
    out[, normal] <- group_patch(out[, normal])
    # patch grouped by treatment
    out[, treatment] <- group_patch(out[, treatment])

    out[, data_cols] <- group_patch(out[, data_cols])  # patch between groups

    # patch normal group by reference
    if (!is.null(ref_data)) {
        out[, data_cols] <- ref_patch(out[, c(id, data_cols)], ref_data)[, -1]
    }

    return(out)
}
