# test_opt.R Copyright 2015 Christian Diener <ch.diener@gmail.com> MIT license.
# See LICENSE for more information.

context("Optimization")

test_that("k optimization works", {
    data(eryth)
    S <- stoichiometry(eryth)
    o <- c(atp = -1, adp = 1)
    i_pyr <- which_reaction(make_irreversible(eryth), S = "adp", P = "pyr")
    expect_true(16 %in% i_pyr$idx)
    opt <- fba(o, S)
    expect_equal(length(opt), ncol(S) + 1)
    expect_true(all(opt <= 1 + sqrt(.Machine$double.eps)))
    expect_error(fba(o, S, v_min = 0.5, v_max = 1))
    p <- opt[-length(opt)] + runif(length(opt)-1)
    cp <- closest(p, S, rep(1, ncol(S)))
    expect_error(closest(p[-1], S, rep(1, ncol(S))))
    expect_true(all(abs(S %*% cp) < sqrt(.Machine$double.eps)))
    expect_true(all(cp > -sqrt(.Machine$double.eps)))
}) 

test_that("objective formulation work", {
    data(eryth)
    S <- stoichiometry(eryth)
    o <- build_objective(eryth[[16]], S)
    expect_equal(length(o), nrow(S))
    expect_true("g6p" %in% names(o))
    expect_error(build_objective(list(NA), S))
    expect_error(build_objective(list(S=1:3, P=1:4, N_S=1, N_P=2), S))
    x <- runif(nrow(S))
    o <- build_objective(x, S)
    expect_equal(x, o) 
})

test_that("parsimonous FBA works", {
    data(eryth)
    S <- stoichiometry(eryth)
    a <- rep(0, ncol(S))
    o <- c(atp = -1, adp = 1)
    opt <- fba(o, S)
    opt_p <- pfba(o, S)
    expect_equal(length(opt_p), ncol(S) + 1)
    expect_true(all(opt_p <= 1))
    expect_true(all(opt_p <= opt + sqrt(.Machine$double.eps)))
    expect_true(abs(opt[16] - opt_p[16]) < sqrt(.Machine$double.eps))
    expect_error(pfba(o, S, v_min = 0.5, v_max = 1))
})
