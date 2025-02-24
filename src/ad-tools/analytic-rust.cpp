#include "../../include/ad-tools/analytic-rust.h"

void init_data_analytic_rust(double **stored_values, int Q) {
    *stored_values = (double *)malloc(Q * NUM_COMPONENTS_STORED_ANALYTIC_RUST * sizeof(double));
}

void free_data_analytic_rust(double **stored_values) {
    if (*stored_values != NULL) {
        free(*stored_values);
        *stored_values = NULL;
    }
}

// Residual Evaluation
void f_analytic_rust(int Q, const double mu, const double lambda, double *dXdx_init, double *dudX, double **stored_values, double *f1) {
    compute_f_analytic(Q, mu, lambda, dXdx_init, dudX, NUM_COMPONENTS_STORED_ANALYTIC_RUST, *stored_values, f1);
}

void df_analytic_rust(int Q, const double mu, const double lambda, double *ddudX, double **stored_values, double *df) {
    compute_df_analytic(Q, mu, lambda, ddudX, NUM_COMPONENTS_STORED_ANALYTIC_RUST, *stored_values, df);
}
