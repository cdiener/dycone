#  test_scraper.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

context("Web scraping")

test_that("SBML scraping works", {
    mod42 = "https://www.ebi.ac.uk/biomodels-main/download?mid=BIOMD0000000042"
    pars = sbml_params(mod42)
    spec = sbml_species(mod42)
    expect_equal(nrow(spec), 15)
    expect_true(!is.null(grep("ATP", spec$name)))    
    expect_equal(nrow(pars), 25)
    expect_true(!is.null(grep("k5f", pars$name)))    
})

test_that("HMDB scraping works", {
    ids = find_hmdb("pyruvate")
    capture.output(concs <- hmdb_concentration(ids[1]))
    expect_true("Pyruvic acid" %in% concs$name)
    expect_more_than(nrow(concs), 1)
})

test_that("patching works", {
    miss = data.frame(id=1:3, n1=NA, n2=c(4,NA,6), d1=NA, d2=c(8,NA,10))
    ref = data.frame(id=2, x=9)
    aim = data.frame(id=1:3, n1=c(4,9,6), n2=c(4,9,6), d1=8:10, d2=8:10)
    full = patch(miss, id=1, 2:3, 4:5, ref_data=ref)
    expect_true(all(aim==full))
})
