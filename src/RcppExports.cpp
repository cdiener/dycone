// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// glpk_fba
NumericVector glpk_fba(IntegerVector ridx, IntegerVector cidx, NumericVector vals, int nrows, int ncols, NumericVector lbs, NumericVector ubs, int objidx);
RcppExport SEXP dycone_glpk_fba(SEXP ridxSEXP, SEXP cidxSEXP, SEXP valsSEXP, SEXP nrowsSEXP, SEXP ncolsSEXP, SEXP lbsSEXP, SEXP ubsSEXP, SEXP objidxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ridx(ridxSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cidx(cidxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vals(valsSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lbs(lbsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ubs(ubsSEXP);
    Rcpp::traits::input_parameter< int >::type objidx(objidxSEXP);
    rcpp_result_gen = Rcpp::wrap(glpk_fba(ridx, cidx, vals, nrows, ncols, lbs, ubs, objidx));
    return rcpp_result_gen;
END_RCPP
}
// glpk_fva
NumericVector glpk_fva(IntegerVector ridx, IntegerVector cidx, NumericVector vals, int nrows, int ncols, NumericVector lbs, NumericVector ubs, int objidx, double objf);
RcppExport SEXP dycone_glpk_fva(SEXP ridxSEXP, SEXP cidxSEXP, SEXP valsSEXP, SEXP nrowsSEXP, SEXP ncolsSEXP, SEXP lbsSEXP, SEXP ubsSEXP, SEXP objidxSEXP, SEXP objfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ridx(ridxSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cidx(cidxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vals(valsSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lbs(lbsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ubs(ubsSEXP);
    Rcpp::traits::input_parameter< int >::type objidx(objidxSEXP);
    Rcpp::traits::input_parameter< double >::type objf(objfSEXP);
    rcpp_result_gen = Rcpp::wrap(glpk_fva(ridx, cidx, vals, nrows, ncols, lbs, ubs, objidx, objf));
    return rcpp_result_gen;
END_RCPP
}
// glpk_pfba
NumericVector glpk_pfba(IntegerVector ridx, IntegerVector cidx, NumericVector vals, int nrows, int ncols, NumericVector lbs, NumericVector ubs, int objidx, double objf);
RcppExport SEXP dycone_glpk_pfba(SEXP ridxSEXP, SEXP cidxSEXP, SEXP valsSEXP, SEXP nrowsSEXP, SEXP ncolsSEXP, SEXP lbsSEXP, SEXP ubsSEXP, SEXP objidxSEXP, SEXP objfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ridx(ridxSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cidx(cidxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vals(valsSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lbs(lbsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ubs(ubsSEXP);
    Rcpp::traits::input_parameter< int >::type objidx(objidxSEXP);
    Rcpp::traits::input_parameter< double >::type objf(objfSEXP);
    rcpp_result_gen = Rcpp::wrap(glpk_pfba(ridx, cidx, vals, nrows, ncols, lbs, ubs, objidx, objf));
    return rcpp_result_gen;
END_RCPP
}