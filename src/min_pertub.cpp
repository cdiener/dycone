/* Copyright 2016 Christian Diener <ch.diener@gmail.com>
   MIT license.
   See LICENSE for more information.
*/

#include <Rcpp.h>
#include <glpk.h>

using namespace Rcpp;

const double ztol = 1.0e-9;


void update_problem(glp_prob* lp, NumericVector cond_terms,
                    NumericVector norm_terms, int neq, double tradeoff) {
    int size = cond_terms.size();
    for(int i=1; i <= size; i++) {
        NumericVector vals = NumericVector::create(cond_terms[i-1],
                                                   norm_terms[i-1]);
        double obj_coef = tradeoff;
        if(is_true(any(is_na(vals)))) {
            vals[0] = 1.0;
            vals[1] = 1.0;
            obj_coef = 1.0 - tradeoff;
        }

        int idx[4] = {0, i, i + size, i + 2*size};
        double coefs[4] = {0.0, 1.0/vals[0], -1.0/vals[1], 1.0};
        glp_set_mat_row(lp, 2*i - 1 + neq, 3, idx, coefs);
        glp_set_row_bnds(lp, 2*i - 1 + neq, GLP_LO, 0.0, 0.0);
        coefs[3] = -1.0;
        glp_set_mat_row(lp, 2*i + neq, 3, idx, coefs);
        glp_set_row_bnds(lp, 2*i + neq, GLP_UP, 0.0, 0.0);
        glp_set_obj_coef(lp, 2*size + i, obj_coef);
    }
}


// [[Rcpp::export]]
SEXP glpk_min_perturb(IntegerVector ridx, IntegerVector cidx,
                   NumericVector vals, int nrows, int ncols,
                   StringVector type, NumericVector row_bounds,
                   NumericVector lbs, NumericVector ubs,
                   NumericMatrix ma_terms, double tradeoff,
                   NumericMatrix perms) {
    int n_vars = ma_terms.nrow();
    glp_term_out(GLP_OFF);
    glp_prob* lp = glp_create_prob();
    glp_set_prob_name(lp, "minimal perturbation");
    glp_add_rows(lp, nrows + 2*n_vars);
    glp_add_cols(lp, ncols);

    double bnd;
    for(int i=1; i <= nrows; i++) {
        bnd = row_bounds[i];
        if (type[i-1] == "equal") {
            glp_set_row_bnds(lp, i, GLP_FX, bnd, bnd);
        } else if (type[i-1] == "larger") {
            glp_set_row_bnds(lp, i, GLP_LO, bnd, bnd);
        } else {
            glp_set_row_bnds(lp, i, GLP_UP, bnd, bnd);
        }
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

    glp_set_obj_dir(lp, GLP_MIN);

    glp_load_matrix(lp, vals.length()-1, ridx.begin(), cidx.begin(),
                    vals.begin());

    glp_smcp params;
    glp_init_smcp(&params);
    params.presolve = GLP_ON;
    params.msg_lev = GLP_MSG_OFF;

    glp_adv_basis(lp, 0);

    NumericMatrix cond_fluxes(n_vars, perms.nrow());
    NumericMatrix norm_fluxes(n_vars, perms.nrow());
    NumericVector obj_vals(perms.nrow());
    Rcout<<"Running permutations";
    for(int p = 0; p < perms.nrow(); p++) {
        Rcout<<".";
        update_problem(lp, ma_terms(_, perms(p, 0) - 1),
                       ma_terms(_, perms(p, 1) - 1), nrows, tradeoff);

        glp_adv_basis(lp, 0);
        int res = glp_simplex(lp, &params);
        if (res != 0) return wrap(res);
        res = glp_get_status(lp);
        if (res != GLP_OPT) return wrap(res);

        for(int i=1; i <= n_vars; i++) {
            cond_fluxes(i - 1, p) = glp_get_col_prim(lp, i);
            norm_fluxes(i - 1, p) = glp_get_col_prim(lp, i + n_vars);
        }
        obj_vals[p] = glp_get_obj_val(lp);
    }
    Rcout<<std::endl;

    glp_delete_prob(lp);

    return List::create(_["disease"] = cond_fluxes,
                        _["normal"] = norm_fluxes,
                        _["obj_values"] = obj_vals);
}
