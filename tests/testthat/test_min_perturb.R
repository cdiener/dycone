# Copyright 2015 Christian Diener <ch.diener@gmail.com>
# MIT license. See LICENSE for more information.

context("minimum perturbation")

test_that("random data gives no significant results", {
    mats <- matrix(rnorm(6 * 68, 10, 1), ncol = 6)
    samples <- factor(rep(c("disease", "normal"), each = 3))
    perturb <- minimal_perturbation(mats, samples, eryth)
    expect_true(all(perturb$k$P.value < 0.05))
    expect_true(all(perturb$k$adj.P.Val > 0.99))
    expect_true(all(perturb$fluxes$P.value < 0.05))
    expect_true(all(perturb$fluxes$adj.P.Val > 0.99))
})
