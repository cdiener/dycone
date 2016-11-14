/* Copyright 2016 Christian Diener <ch.diener@gmail.com>
   MIT license.
   See LICENSE for more information.
*/

#include <Rcpp.h>
#include <glpk.h>

using namespace Rcpp;

const double ztol = 1.0e-9;

// [[Rcpp::export]]
NumericVector glpk_fva(IntegerVector ridx, IntegerVector cidx,
                   NumericVector vals, int nrows, int ncols,
                   NumericVector lbs, NumericVector ubs,
                   int objidx, double objf) {
    glp_term_out(GLP_OFF);
    glp_prob* lp = glp_create_prob();
    glp_set_prob_name(lp, "FBA problem");
    glp_add_rows(lp, nrows);
    glp_add_cols(lp, ncols);

    for(int i=1; i <= nrows; i++) {
        glp_set_row_bnds(lp, i, GLP_FX, 0.0, 0.0);
    }

    int eq = 0;
    for(int i=1; i <= ncols; i++) {
        eq = (ubs[i] - lbs[i]) < ztol;
        if (eq) {
            glp_set_col_bnds(lp, i, GLP_FX, lbs[i], ubs[i]);
        } else {
            glp_set_col_bnds(lp, i, GLP_DB, lbs[i], ubs[i]);
        }
    }

    glp_set_col_bnds(lp, objidx, GLP_LO, 0.0, 0.0);
    glp_set_obj_coef(lp, objidx, 1.0);
    glp_set_obj_dir(lp, GLP_MAX);

    glp_load_matrix(lp, vals.length()-1, ridx.begin(), cidx.begin(),
                    vals.begin());

    glp_smcp params;
    glp_init_smcp(&params);
    params.presolve = GLP_ON;
    params.msg_lev = GLP_MSG_OFF;

    glp_adv_basis(lp, 0);
    int res = glp_simplex(lp, &params);
    if (res != 0) return res;
    res = glp_get_status(lp);
    if (res != GLP_OPT) return res;
    double obj_sol = glp_get_obj_val(lp);

    NumericMatrix ranges(ncols, 2);
    params.presolve = GLP_OFF;
    glp_set_col_bnds(lp, objidx, GLP_LO, objf*obj_sol, objf*obj_sol);
    int dir[2] = {GLP_MIN, GLP_MAX};
    glp_set_obj_coef(lp, objidx, 0.0);
    for(int col=0; col<2; col++) {
        for(int i=1; i <= ncols; i++) {
            if (i == objidx) {
                ranges(i-1, col) = objf*obj_sol;
                continue;
            }
            glp_set_obj_coef(lp, i, 1.0);
            glp_set_obj_dir(lp, dir[col]);
            res = glp_simplex(lp, &params);
            if (res != 0) return res;
            res = glp_get_status(lp);
            if (res != GLP_OPT) return res;
            ranges(i-1, col) = glp_get_obj_val(lp);
            glp_set_obj_coef(lp, i, 0.0);
        }
    }

    glp_delete_prob(lp);

    return ranges;
}
