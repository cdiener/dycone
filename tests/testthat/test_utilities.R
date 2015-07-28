#  utility_tests.R
#  
#  Copyright 2015 Christian Diener <ch.diener@gmail.com>
#  
#  MIT license. See LICENSE for more information.

context("Model utilities")

test_that("we can convert extra information nicely", {
    expect_equal(str_conv("1"), 1)
    expect_equal(str_conv("1,  2"), 1:2)
    expect_equal(str_conv("a, 1"), c("a","1"))
    expect_equal(str_conv(" abc "), "abc")
    expect_equal(str_conv("abc,gef"), c("abc", "gef"))
})

test_that("order_by works", {
    expect_equal(order_by(1:3, 3:1), 3:1)
})

test_that("the caching operator works", {
    file.remove("test-42.Rd")
    { x <- 1:10; y <- x*x } %c% "test-42.Rd"
    {x <- 42; y <- "bla"} %c% "test-42.Rd"
    load("test-42.Rd")
    file.remove("test-42.Rd")
    expect_equal(x,1:10)
    expect_equal(y, (1:10)^2)
})

test_that("we can calculate mass-action terms", {
    expect_equal(mass_action(c(-2,3), 2:3), 4)
    expect_equal(mass_action(c(-1,3,-2), 1:3), 9)
    expect_equal(mass_action(c(2,3), 2:3), 1)
})

test_that("we can calculate derivatives of mass-action terms", {
    expect_equal(partial_deriv(1, c(-2,3), 2:3), 2)
    expect_equal(partial_deriv(3, c(-1,3,-2), 1:3), 3)
    expect_equal(partial_deriv(1, c(2,3), 2:3), 1)
})

test_that("we can parse reactions with no substrates or products", {
	expect_true( is.na(get_reaction_elems("-> B")$S) )
	expect_true( is.na(get_reaction_elems("A ->")$P) )
	expect_equal(get_reaction_elems("-> B")$N_S, 1)
	expect_equal(get_reaction_elems("A ->")$N_P, 1)
})

test_that("we can identify reversibility", {
	expect_true( get_reaction_elems("A <-> B")$rev )
	expect_true( get_reaction_elems("A <=> B")$rev )
	expect_true( get_reaction_elems("A < - > B")$rev )
	expect_false( get_reaction_elems("A -> B")$rev )
	expect_false( get_reaction_elems("A => B")$rev )
})

test_that("we can get info from sample model", {
	data(eryth)
	s = get_species(eryth)
	n_s = length(s)
	expect_true( all(grepl("\\d:", format(eryth))) )
	expect_true( all(dim(get_stochiometry(eryth)) == c(n_s-1, 68)) )
	expect_true( all(dim(get_stochiometry(eryth, reversible=TRUE)) == c(n_s-1, 45)) )
	expect_true( all(c("atp", "23dpg", "none") %in% get_species(eryth)) )
})
