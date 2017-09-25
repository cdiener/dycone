# Copyright 2017 Christian Diener <ch.diener@gmail.com>
#
# MIT license. See LICENSE for more information.

#' Assemble the LP problem for minimal perturbation.
#'
#' @param ma_terms A matrix or data frame containing the metabolic terms in the
#'  columns.
#' @param samples A factor or character string with either "normal" and
#'  "disease" entries.
#' @param S The stoichiometric matrix or a list containing one matrix for each
#'  sample.
#' @param bounds The bounds for the fluxes. Either a single value for the max
#'  bound or a list each containing a vector specifying the bounds for the
#'  specific sample.
#' @param tradeoff Double in [0, 1] specifiying the equlibrium between
#'  minimizing differences in k changes versus flux changes.
#' @param growth Constraints fro growth rate. If not NA must be a named vector
#'  with entries "idx", "min", "max" denoting the index of the biomass reaction
#'  and its minimum and maximum flux respectively.
#' @return A list with the following components:
#'  \describe{
#'  \item{coefficients}{The coefficients for the equalities/inequalities in
#'      the LP. Each row denotes one equation.}
#'  \item{type}{The type of the (in)equalities. A character vector with as many
#'      entries as rows in `coefficients` being either "equal", "less" or
#'      "larger"}
#'  \item{bounds}{A matrix with two columns containing the upper and lower
#'      variable (flux) bounds.}
#'  \item{row_bounds}{Vector containing the right sides to the (in)equalities.}
#'  \item{obj_coefs}{The objective coefficients. Containing one number for each
#'      variable.}
#'  }
#' @examples
#'  NULL
#'
#' @export
#' @importFrom Matrix Matrix bdiag
minimal_perturbation_lp <- function(ma_terms, samples, S, v_min = 0,
                                    v_max = 1000, min_obj = 100) {
    if (!is.list(S)) {
        new_S <- vector(mode = "list", length = 2)
        for (i in 1:2) {
            new_S[[i]] <- S
        }
        S <- new_S
    }
    n_var <- ncol(S[[1]])
    coefficients <- bdiag(S)
    coefficients <- cbind(coefficients, Matrix(0, ncol = n_var,
                          nrow = nrow(coefficients)))
    type <- rep("equal", nrow(coefficients))
    row_bounds <- rep(0, nrow(coefficients))

    if (is.list(v_min)) {
        v_min <- do.call(c, v_min)
    } else {
        v_min <- rep(v_min, ncol(coefficients))
    }
    if (is.list(v_max)) {
        v_max <- do.call(c, v_min)
    } else {
        v_max <- rep(v_max, ncol(coefficients))
    }
    if (length(min_obj) == 2) {
        idx <- (0:1) * n_vars + min_obj[1]
        v_min[idx] <- min_obj[2]
        v_max[idx] <- pmax(v_max[idx], min_obj[2])
    } else {
        coefs <- matrix(0, nrow = 2, ncol = 3 * n_var)
        for (i in 1:2) {
            coefs[i, ((i - 1) * n_var + 1):(i * n_var)] <- 1
        }
        coefficients <- rbind(coefficients, coefs)
        type <- append(type, rep("larger", 2))
        row_bounds <- append(row_bounds, rep(min_obj, 2))
    }

    return(list(
        coefficients = coefficients,
        type = type,
        bounds = matrix(c(v_min, v_max), ncol = 2),
        row_bounds = row_bounds
    ))
}

#' Assemble the LP problem for minimal perturbation.
#'
#' @param ma_terms A matrix or data frame containing the metabolic terms in the
#'  columns.
#' @param samples A factor or character string with either "normal" and
#'  "disease" entries.
#' @param reacts The reactions for the respective model.
#' @param bounds The bounds for the fluxes. Either a single value for the max
#'  bound or a list each containing a vector specifying the bounds for the
#'  specific sample.
#' @param tradeoff Double in [0, 1] specifiying the equlibrium between
#'  minimizing differences in k changes versus flux changes.
#' @param growth Constraints fro growth rate. If not NA must be a named vector
#'  with entries "idx", "min", "max" denoting the index of the biomass reaction
#'  and its minimum and maximum flux respectively.
#' @return A list with the following components:
#'  \describe{
#'  \item{fluxes}{The alterations in fluxes.}
#'  \item{k}{The alterations in kinetic constants kcat.}
#'  }
#' @examples
#'  NULL
#'
#' @export
#' @importFrom Matrix Matrix bdiag
minimal_perturbation <- function(ma_terms, samples, reacts,
                                 permutations = 1000, v_min = 0,
                                 v_max = 1000, tradeoff = 0.95,
                                 min_obj = 100) {
    S <- stoichiometry(reacts, sparse = TRUE)
    prob <- minimal_perturbation_lp(ma_terms, samples, S, v_min, v_max,
                                    min_obj)
    sp <- summary(prob$coefficients)
    idxs <- list(which(samples == levels(samples)[1]),
                 which(samples == levels(samples)[2]))
    if (permutations >= prod(table(samples))) {
        perms <- do.call(expand.grid, idxs)
    } else {
        perms <- do.call(cbind, lapply(idxs, sample,
                         size = permutations, replace = TRUE))
    }
    perms <- unname(as.matrix(perms))

    sol <- glpk_min_perturb(
        gf(sp[,1]), gf(sp[,2]), gf(sp[,3]), nrow(prob$coefficients),
        ncol(prob$coefficients), prob$type, gf(prob$row_bounds),
        gf(prob$bounds[,1]), gf(prob$bounds[,2]), ma_terms,
        0.95, perms)

    if (!is.list(sol)) {
        stop(paste0("GLPK failed with error code ", sol, "."))
    }

    samples <- factor(rep(c("normal", "disease"),
                      times = c(ncol(sol$normal), ncol(sol$disease))))
    fluxes <- cbind(sol$normal, sol$disease)
    fluxes[fluxes < 0] <- 0
    ks <- fluxes / ma_terms[, c(perms[, 2], perms[, 1])]

    design <- model.matrix(~ 0 + samples)
    colnames(design) <- levels(samples)
    contrasts <- makeContrasts(disease - normal, levels = design)
    flux_fit <- contrasts.fit(lmFit(fluxes, design), contrasts)
    flux_fit <- eBayes(flux_fit)
    k_fit <- contrasts.fit(lmFit(ks, design), contrasts)
    k_fit <- eBayes(k_fit)

    reacts <- make_irreversible(reacts)
    ids <- sapply(reacts, `[[`, "abbreviation")
    rstr <- sapply(reacts, to_string, name = FALSE)

    flux_stats <- topTable(flux_fit, number = nrow(fluxes), confint = TRUE,
                           sort.by = "none", adjust.method = "fdr")
    flux_stats <- cbind(data.frame(id = ids, reaction = rstr), flux_stats)
    k_stats <- topTable(k_fit, number = nrow(ks), confint = TRUE,
                           sort.by = "none", adjust.method = "fdr")
    k_stats <- cbind(data.frame(id = ids, reaction = rstr)[complete, ],
                     k_stats)

    return(list(fluxes = flux_stats, k = k_stats,
                raw = fluxes, obj_value = sol$obj_value))
}
