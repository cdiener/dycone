# test_opt.R Copyright 2015 Christian Diener <ch.diener@gmail.com> MIT license.
# See LICENSE for more information.

context("Optimization")

test_that("k optimization works", {
    data(eryth)
    S <- stochiometry(eryth)
    a <- rep(0, ncol(S))
    i_pyr <- which_reaction(make_irreversible(eryth), S = "adp", P = "pyr")
    expect_true(16 %in% i_pyr$idx)
    a[16] <- 1
    opt <- dba(a, S, rep(1, ncol(S)), lower = 0, upper = 1)
    expect_equal(length(opt), ncol(S))
    expect_true(all(opt <= 1))
    expect_error(dba(a, S, rep(1, ncol(S)), lower = 0.5, upper = 1))
    p <- opt + runif(length(opt))
    cp <- closest(p, S, rep(1, ncol(S)))
    expect_error(closest(p[-1], S, rep(1, ncol(S))))
    expect_true(all(abs(S %*% cp) < sqrt(.Machine$double.eps)))
    expect_true(all(cp > -sqrt(.Machine$double.eps)))
}) 

test_that("parsimonous FBA works", {
    data(eryth)
    S <- stochiometry(eryth)
    a <- rep(0, ncol(S))
    a[16] <- 1
    i_pyr <- which_reaction(make_irreversible(eryth), S = "adp", P = "pyr")
    opt <- dba(a, S, rep(1, ncol(S)), lower = 0, upper = 1)
    opt_p <- pba(a, S, rep(1, ncol(S)), lower = 0, upper = 1)
    expect_equal(length(opt_p), ncol(S))
    expect_true(all(opt_p <= 1))
    expect_true(all(opt_p <= opt + sqrt(.Machine$double.eps)))
    expect_true(abs(opt[16] - opt_p[16]) < sqrt(.Machine$double.eps))
    expect_error(pba(a, S, rep(1, ncol(S)), lower = 0.5, upper = 1))
})
