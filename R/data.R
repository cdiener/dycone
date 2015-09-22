# data.R Copyright 2015 Christian Diener <ch.diener@gmail.com> MIT license. See
# LICENSE for more information.

#' dycone: Differential k-cone analysis for metabolome data.
#'
#' @docType package
#' @name dycone
#' @importFrom magrittr '%>%'
NULL

### Data documentation

#' Metabolic network of human red blood cell.
#'
#' Contains the reactions of the human red blood cell (erythrocite) metabolic 
#' model together with some randomly sampled rates k.
#'
#' @return Assigns the \code{eryth} reactions list to the namespace.
#' @format A list of reactions.
#' @source \url{http://journals.plos.org/plosone/article/asset?unique&id=
#' info:doi/10.1371/journal.pone.0004967.s003}
"eryth" 
