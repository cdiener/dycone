# test_kcone.R Copyright 2015 Christian Diener <ch.diener@gmail.com> MIT license.
# See LICENSE for more information.

context("kcone analysis")

test_that("cones can be calculated and analyzed", {
    S <- t(c(1, -1))
    rownames(S) <- "A"
    V <- polytope_basis(S)
    expect_equal(as.vector(V), c(1, 1))
    expect_equal(as.vector(kcone_null(S, c(1, 1))), c(sqrt(2), sqrt(2))/2)
    mats <- runif(2, 0, 10)
    K <- kcone(V, mats)
    expect_equal(dim(K), dim(V))
    expect_equal(K[2]/K[1], mats[1]/mats[2])
    K <- kcone(V, mats, normalize = T)
    expect_equal(K[2]/K[1], mats[1]/mats[2])
    expect_equal(sum(K * K), 1)
    expect_true(inside(rev(mats * 0.42), S, mats))
})

test_that("stability analysis works", {
    S <- matrix(c(1, 0, 0, 1, -1, 0, 0, -1), nrow = 2)
    rownames(S) <- c("A", "B")
    V <- polytope_basis(S)
    stab <- stability_analysis(V, S, c(A = 1, B = 1))
    expect_true(all(stab$what == "stable"))
    expect_equal(stab$ev1, c(0, 0))
    expect_equal(stab$ev2, -0.5 * rep(sqrt(2), 2))
})

test_that("helper functions work", {
    x <- runif(10)
    x <- x/enorm(x)
    expect_equal(sum(x * x), 1)
    B <- matrix(c(0, 1, 1, 1), ncol = 2)
    class(B) <- "basis"
    expect_equal(occupation(B), c(0.5, 1))
    x <- 1:0
    y <- c(0, 2)
    expect_equal(angle(x, y), pi/2)
    expect_equal(scaling(y, 0:1), 2)
})

test_that("eigendynamics work", {
    B <- matrix(runif(10000), nrow = 10)
    class(B) <- append(class(B), "basis")
    expect_true(all(eigendynamics(B) >= 0))
    expect_equal(dim(eigendynamics(B, 3)), c(nrow(B), 3))
})

test_that("hypothesis generation works", {
    data(eryth)
    n_r <- length(make_irreversible(eryth))
    mats <- matrix(runif(6 * n_r), ncol = 6)
    h <- hyp(mats[, 1:3], mats[, 4:6], eryth)
    expect_equal(nrow(h), n_r)
    expect_true(all(h$pval <= 1))
    h_opt <- hyp(mats[, 1:3], mats[, 4:6], eryth, type = "optimization", v_opt = runif(n_r), 
        full = T)
    expect_equal(length(h_opt), 7)
    expect_equal(length(h_opt$obj_normal), 3)
    expect_equal(length(h_opt$obj_disease), 3)
    expect_equal(ncol(h_opt$k_normal), 3)
    expect_equal(ncol(h_opt$k_disease), 3)
    expect_equal(nrow(h_opt$hyp), n_r)
})

test_that("basis plotting works", {
    B <- matrix(runif(10000), nrow = 10)
    class(B) <- append(class(B), "basis")
    png("test-plot-42.png")
    out <- capture.output(plot_basis(B, n_cl = 100))
    dev.off()
    expect_true(file.exists("test-plot-42.png"))
    expect_match(out, "in-cluster distance", all = F)
    file.remove("test-plot-42.png")
    png("test-plot-42.png")
    out <- capture.output(plot_red(list(B), n_cl = 100))
    dev.off()
    expect_true(file.exists("test-plot-42.png"))
    expect_match(out, "Information", all = F)
    expect_match(out, "in-cluster", all = F)
    file.remove("test-plot-42.png")
}) 
