# test_sbml.R Copyright 2015 Christian Diener <ch.diener@gmail.com> MIT
# license. See LICENSE for more information.

context("SBML")

test_that("SBML scraping works", {
    mod42 <- "http://www.ebi.ac.uk/biomodels-main/download?mid=BIOMD0000000042"
    pars <- sbml_params(mod42)
    spec <- sbml_species(mod42)
    expect_output(reacts <- sbml_reactions(mod42), "=")
    expect_equal(nrow(spec), 15)
    expect_true(!is.null(grep("ATP", spec$name)))
    expect_equal(nrow(pars), 25)
    expect_true(!is.null(grep("k5f", pars$name)))
    expect_equal(length(reacts), 25)
})

test_that("fbc bounds can be read", {
    cobra_minifbc <- "https://raw.githubusercontent.com/opencobra/cobrapy/master/cobra/test/data/mini_fbc2.xml"
    out <- capture.output(mod <- read_sbml(cobra_minifbc))
    expect_true(length(grep("with FBC", out)) > 0)
    expect_true(length(grep("reactions mapped", out)) > 0)
    mapped <- sapply(mod$reactions, function(x) 
        all(is.numeric(c(x$lower, x$upper))))
    expect_true(all(mapped))
})

test_that("SBML writing works", {
    data(eryth)
    s <- species(eryth)
    s <- data.frame(id=s, kegg.compound="CXXXXX", hmdb="HMDBXXXXX")
    am <- c(abbreviation="dummy")
    out <- tempfile("model", fileext=".xml")
    expect_silent(write_sbml(eryth, s, annmap=am, out=out))
    read_m <- capture.output(r <- read_sbml(out))
    expect_equal(nrow(r$species), nrow(s))
    expect_equal(length(r$reactions), length(eryth))
    expect_true("kegg.compound" %in% names(r$species))
    expect_true("hmdb" %in% names(r$species))
    expect_true("dummy" %in% names(r$reactions[[1]]))
})
