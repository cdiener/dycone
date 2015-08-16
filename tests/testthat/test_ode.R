# test_ode.R Copyright 2015 Christian Diener <ch.diener@gmail.com> MIT license.
# See LICENSE for more information.

context("ODE solving")

test_that("timecourse works", {
    data(eryth)
    S <- stochiometry(eryth)
    x0 <- runif(nrow(S), 0, 1)
    names(x0) <- rownames(S)
    tc <- timecourse(x0, seq(0, 10, 0.1), runif(ncol(S), 0, 1), S)
    expect_equal(tc[, 1], seq(0, 10, 0.1))
    expect_true(all(names(x0) %in% colnames(tc)))
}) 
