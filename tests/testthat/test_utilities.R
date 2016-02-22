# test_utilities.R Copyright 2015 Christian Diener <ch.diener@gmail.com> MIT
# license. See LICENSE for more information.

context("Model utilities")

test_that("we can convert extra information nicely", {
    expect_equal(str_conv("1"), 1)
    expect_equal(str_conv("1,  2"), 1:2)
    expect_equal(str_conv("a, 1"), c("a", "1"))
    expect_equal(str_conv(" abc "), "abc")
    expect_equal(str_conv("abc,gef"), c("abc", "gef"))
})

test_that("orderby works", {
    expect_equal(orderby(1:3, 3:1), 3:1)
})

test_that("the caching operator works", {

    f <- tempfile("cache", fileext=".Rd")
    {
        x <- 1:10
        y <- x * x
    } %c% f
    {
        x <- 42
        y <- "bla"
    } %c% f
    load(f)
    expect_equal(x, 1:10)
    expect_equal(y, (1:10)^2)
})

test_that("we can calculate mass-action terms", {
    expect_equal(mass_action(c(-2, 3), 2:3), 4)
    expect_equal(mass_action(c(-1, 3, -2), 1:3), 9)
    expect_equal(mass_action(c(2, 3), 2:3), 1)
})

test_that("we can calculate all mass-action terms", {
    co <- c(x = 1, y = 2)
    S <- matrix(c(-1, 1, 3, -2), ncol = 2, byrow = T)
    rownames(S) <- c("x", "y")
    expect_equal(ma_terms(S, co), c(1, 4))
    expect_error(ma_terms(S, 1:2))
    
    rownames(S) <- c("y", "x")
    expect_equal(ma_terms(S, co), 2:1)
    
    co <- data.frame(name = c("x", "y"), co, co)
    expect_equal(rowSums(ma_terms(S, co)), c(4, 2))
    expect_error(ma_terms(S, co[-1]))
})

test_that("we can cancalculate a Jacobian", {
    S <- matrix(c(-1, 1, 3, -2), ncol = 2, byrow = TRUE)
    concs <- c(2, 4)
    rownames(S) <- names(concs) <- c("A", "B")
    J <- matrix(c(1, 0, 0, 8), ncol = 2, byrow = TRUE)
    expect_equal(jacobian(S, concs), J)
})

test_that("we can calculate derivatives of mass-action terms", {
    expect_equal(deriv_ma(1, c(-2, 3), 2:3), 4)
    expect_equal(deriv_ma(3, c(-1, 3, -2), 1:3), 6)
    expect_equal(deriv_ma(1, c(2, 3), 2:3), 0)
})

test_that("we can parse reactions with no substrates or products", {
    expect_true(is.na(get_reaction_elems("-> B")$S))
    expect_true(is.na(get_reaction_elems("A ->")$P))
    expect_equal(get_reaction_elems("-> B")$N_S, 1)
    expect_equal(get_reaction_elems("A ->")$N_P, 1)
})

test_that("we can identify reversibility", {
    expect_true(get_reaction_elems("A <-> B")$rev)
    expect_true(get_reaction_elems("A <=> B")$rev)
    expect_true(get_reaction_elems("A < - > B")$rev)
    expect_false(get_reaction_elems("A -> B")$rev)
    expect_false(get_reaction_elems("A => B")$rev)
})

test_that("we can read and format reactions from a file", {
    r_str <- "reaction,abbreviation,numbers\nA -> B,blub,\"1,2,3\"\nB <=>, bla, 3"
    r <- read_reactions(textConnection(r_str))
    expect_equal(length(r), 2)
    expect_equal(r[[1]]$numbers, 1:3)
    expect_equal(r[[2]]$abbreviation, "bla")
    expect_false(r[[1]]$rev)
    expect_true(r[[2]]$rev)
    expect_equal(species(r), c("A", "B"))
    expect_equal(grep("1*A", format(r)), 1)
    expect_equal(grep("Model has", capture.output(print(r))), 1)
    
    r_str <- "reaction,abbreviation,numbers\nA -Z B,blub,\"1,2,3\"\nB <=>, bla, 3"
    expect_error(invisible(read_reactions(textConnection(r_str))), "arrows")
})

test_that("we can get the stoichiometry from reactions", {
    r_str <- "reaction,abbreviation,numbers\nA -> B,blub,\"1,2,3\"\nB <=>, bla, 3"
    r <- read_reactions(textConnection(r_str))
    S1 <- matrix(c(-1, 1, 0, -1, 0, 1), ncol = 3)
    S2 <- matrix(c(-1, 1, 0, -1), ncol = 2)
    S1s <- Matrix::Matrix(c(-1, 1, 0, -1, 0, 1), sparse = TRUE, ncol = 3)
    rownames(S1s) <- c("A", "B")
    rownames(S1) <- rownames(S2) <- c("A", "B")
    expect_equal(stoichiometry(r), S1)
    expect_equal(stoichiometry(r, sparse=TRUE), S1s)
    expect_equal(stoichiometry(r, reversible = TRUE), S2)
})

test_that("we can convert to irreversible", {
    r_str <- "reaction,abbreviation,numbers\nA -> B,blub,\"1,2,3\"\nB <=>, bla, 3"
    r <- read_reactions(textConnection(r_str))
    expect_equal(length(make_irreversible(r)), 3)
})

test_that("we can get info from sample model", {
    data(eryth)
    s <- species(eryth)
    n_s <- length(s)
    expect_true(all(grepl("\\d:", format(eryth))))
    expect_true(all(dim(stoichiometry(eryth)) == c(n_s, 68)))
    expect_true(all(dim(stoichiometry(eryth, reversible = TRUE)) == c(n_s, 45)))
    expect_true(all(c("atp", "23dpg") %in% species(eryth)))
})

test_that("graph conversions", {
    r <- list(a = "throws", b = "error")
    expect_error(as.graph(r))
    data(eryth)
    expect_true(class(as.graph(eryth)) %in% c("matrix", "igraph"))
    png("test-plot-42.png")
    plot(eryth)
    dev.off()
    expect_true(file.exists("test-plot-42.png"))
    file.remove("test-plot-42.png")
})

test_that("additional helpers work", {
    r_str <- "reaction,abbreviation,numbers\nA -> B,blub,\"1,2,3\"\nB <=>, bla, 3"
    r <- read_reactions(textConnection(r_str))
    expect_equal(r_order(r), c(1, 1))
    expect_equal(r_order(make_irreversible(r)), c(1, 1, 0))
    expect_false(any(constant_flux(r)))
    expect_true(any(constant_flux(make_irreversible(r))))
    expect_true("reactions" %in% class(as.reactions(stoichiometry(r))))
    expect_more_than(length(methods(class = "reactions")), 1)
    expect_equal(rp(r, "numbers")[, 2], c(1, 2, 3, 3))
}) 
