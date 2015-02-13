#  utility_tests.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

context("Model utilities")

test_that("can parse reactions with no substrates or products", {
	expect_true( is.na(get_reaction_elems("-> B")$S) )
	expect_true( is.na(get_reaction_elems("A ->")$P) )
	expect_equal(get_reaction_elems("-> B")$N_S, 1)
	expect_equal(get_reaction_elems("A ->")$N_P, 1)
})

test_that("can identify reversibility", {
	expect_true( get_reaction_elems("A <-> B")$rev )
	expect_true( get_reaction_elems("A <=> B")$rev )
	expect_true( get_reaction_elems("A < - > B")$rev )
	expect_false( get_reaction_elems("A -> B")$rev )
	expect_false( get_reaction_elems("A => B")$rev )
})

test_that("can get info from sample model", {
	data(eryth)
	s = get_species(eryth)
	n_s = length(s)
	expect_true( all(grepl("\\d:", format(eryth))) )
	expect_true( all(dim(get_stochiometry(eryth)) == c(n_s-1, 68)) )
	expect_true( all(dim(get_stochiometry(eryth, reversible=TRUE)) == c(n_s-1, 45)) )
	expect_true( all(c("atp", "23dpg", "none") %in% get_species(eryth)) )
})
