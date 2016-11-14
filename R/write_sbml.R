#  write_sbml.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

HEADER <- '<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2" level="3" sboTerm="SBO:0000624" version="1" xmlns="http://www.sbml.org/sbml/level3/version1/core" fbc:required="false">
  <model fbc:strict="true" id="%s">
    <listOfUnitDefinitions>
      <unitDefinition id="mmol_per_gDW_per_hr">
        <listOfUnits>
          <unit exponent="1" kind="mole" multiplier="1" scale="-3"/>
          <unit exponent="-1" kind="gram" multiplier="1" scale="0"/>
          <unit exponent="-1" kind="second" multiplier="3600" scale="0"/>
        </listOfUnits>
      </unitDefinition>
    </listOfUnitDefinitions>
    <listOfCompartments>
      <compartment constant="true" id="cell" name="cellular interior"/>
    </listOfCompartments>'

OBJECTIVE <- '    <fbc:listOfObjectives fbc:activeObjective="obj">
      <fbc:objective fbc:id="obj" fbc:type="maximize">
        <fbc:listOfFluxObjectives>
          <fbc:fluxObjective fbc:reaction="%s" fbc:coefficient="1"/>
        </fbc:listOfFluxObjectives>
      </fbc:objective>
    </fbc:listOfObjectives>'

SPECIES <- '      <species boundaryCondition="false" constant="false" hasOnlySubstanceUnits="false" id="%s" name="%s" metaid="%s" initialConcentration="0" compartment="cell">'

SPEC_REF <- '          <speciesReference constant="true" species="%s" stoichiometry="%f"/>'

RDF_OPEN <- '        <sbml:annotation xmlns:sbml="http://www.sbml.org/sbml/level3/version1/core">
          <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
            <rdf:Description rdf:about="#%s">
              <bqbiol:is xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
                <rdf:Bag>'

RDF_CLOSE <- '                </rdf:Bag>
              </bqbiol:is>
            </rdf:Description>
          </rdf:RDF>
        </sbml:annotation>'

RDF_LI <- '                  <rdf:li rdf:resource="http://identifiers.org/%s/%s"/>'

REACT <- '      <reaction fast="false" id="%s" reversible="%s" fbc:lowerFluxBound="%s" fbc:upperFluxBound="%s">'

DPARAM <- '      <parameter constant="true" id="%s" sboTerm="SBO:0000626" units="mmol_per_gDW_per_hr" value="%f"/>'

PARAM <- '      <parameter constant="true" id="%s" sboTerm="SBO:0000625" units="mmol_per_gDW_per_hr" value="%f"/>'

sid <- function(id) {
    paste0("M_", gsub("[^\\w\\d_]", "_", id, perl = T))
}

write_param <- function(r, i, con) {
    if (!is.null(r$lower))
        write(sprintf(PARAM, paste0("r", i, "_lower"), r$lower), file = con)
    if (!is.null(r$upper))
        write(sprintf(PARAM, paste0("r", i, "_upper"), r$upper), file = con)
}

write_species <- function(spec, ann, con) {
    s <- sid(spec)
    write(sprintf(SPECIES, s, s, s), file = con)
    if (any(!is.na(ann))) {
        write(sprintf(RDF_OPEN, s), file = con)
        for (i in 1:length(ann)) {
            if (!is.na(ann[i])) {
                x <- unique(as.character(str_conv(ann[i])))
                nx <- names(ann)[i]
                sapply(x, function(id) write(sprintf(RDF_LI, nx, id),
                                             file = con))
            }
        }
        write(RDF_CLOSE, file = con)
    }
    write("      </species>", file = con)
}

write_reaction <- function(r, i, db, annmap, con) {
    rev_str <- c("false", "true")[r$rev + 1]
    name <- paste0("r", i)
    lower <- paste0(name, "_lower")
    upper <- paste0(name, "_upper")
    if (is.null(r$lower)) {
        if (r$rev) lower <- "default_lower" else lower <- "zero_bound"
    }
    if (is.null(r$upper)) upper <- "default_upper"
    write(sprintf(REACT, name, rev_str, lower, upper), file = con)

    # Write annotations
    val <- names(r)[names(r) %in% names(annmap)]
    if (length(val) > 0 & !all(is.na(r[val]))) {
        write(sprintf(RDF_OPEN, name), file = con)
        for (v in val) {
             if (!all(is.na(r[[v]]))) {
                x <- unique(as.character(r[[v]]))
                nx <- annmap[v]
                sapply(x, function(id) write(sprintf(RDF_LI, nx, id),
                                             file = con))
            }
        }
        write(RDF_CLOSE, file = con)
    }

    # Reactants and substrates
    if (!is.na(r$S)[1]) {
        write("        <listOfReactants>", file = con)
        idx <- 1:length(r$S)
        sapply(idx, function(i)
            write(sprintf(SPEC_REF, sid(r$S[i]), r$N_S[i]), file = con))
        write("        </listOfReactants>", file = con)
    }

    # Products
    if (!is.na(r$P)[1]) {
        idx <- 1:length(r$P)
        write("        <listOfProducts>", file = con)
        sapply(idx, function(i)
            write(sprintf(SPEC_REF, sid(r$P[i]), r$N_P[i]), file = con))
        write("        </listOfProducts>", file = con)
    }

    write("      </reaction>", file = con)
}

#" Writes a SBML model.
#"
#" This function writes a reaction list as a SBML file. It accepts various levels
#" of additional information to help annotate the model. \code{write_sbml} will
#" only write SBML Level 3 Version 1 with FBC Version 2.
#"
#" @seealso \code{\link{read_sbml}}, to read SBML models.
#" @param reacts A reaction list.
#" @param spec Additional annotations for the species. Must be a data frame where
#"  the first column denotes the name/id and additional columns are valid id
#"  types from http://identifiers.org.
#" @param obj A positive integer denoting the index of the reaction taken to be
#"  the optimization target. If NA no FBA objective will be added.
#" @param default_bounds numeric vector of length 2. The default bounds used if
#"  no explicit lower and upper bounds are given for the reaction. If the reaction
#"  is irreversible the lower bounds will be 0.
#" @param annmap An annotation map. A named character vector where the names denotes
#"  reactions annotations in the reation list and the entries denote the respective
#"  id on http://identifiers.org.
#" @param out A string denoting the name of the written SBML file.
#" @return Nothing.
#" @examples
#" data(eryth)
#" write_sbml(eryth, out="eryth.xml")
#"
#" @export
write_sbml <- function(reacts, spec=NA, obj=NA, default_bounds=c(-1000, 1000),
    annmap=NA, out="dycone_model.xml") {
    if (is.null(ncol(spec))) {
        spec <- species(reacts)
        spec <- data.frame(id = spec, ann = NA)
    } else if (!all(species(reacts) %in% spec[, 1]))
        stop(paste0("There are species in the reactions which are not",
                    "contained in 'spec'!"))

    con <- file(out, open = "w")

    # write the SBML header
    write(sprintf(HEADER, sub("\\..+", "", basename(out))), file = con)

    # write the FBA objective if given
    if (!is.na(obj)) write(sprintf(OBJECTIVE, paste0("r", obj)), file = con)

    # write the list of parameters
    write("    <listOfParameters>", file = con)
    write(sprintf(DPARAM, "zero_bound", 0), file = con)
    write(sprintf(DPARAM, "default_lower", default_bounds[1]), file = con)
    write(sprintf(DPARAM, "default_upper", default_bounds[2]), file = con)
    d <- sapply(1:length(reacts), function(i) write_param(reacts[[i]], i, con))
    write("    </listOfParameters>", file = con)

    # write the species list
    write("    <listOfSpecies>", file = con)
    dummy <- apply(spec, 1, function(s) write_species(s[1], s[-1], con))
    write("    </listOfSpecies>", file = con)

    # write the reaction list
    write("    <listOfReactions>", file = con)
    dummy <- sapply(1:length(reacts), function(i)
        write_reaction(reacts[[i]], i, default_bounds, annmap, con))
    write("    </listOfReactions>", file = con)

    write("  </model>\n</sbml>", file = con)
    close(con)
}
