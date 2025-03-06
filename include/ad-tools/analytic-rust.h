#ifndef ANALYTIC_RUST_H
#define ANALYTIC_RUST_H

#include "../utils.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NUM_COMPONENTS_STORED_ANALYTIC_RUST 15

void init_data_analytic_rust(double **stored_values, int Q, int *num_comp);
void free_data_analytic_rust(double **stored_values);
void f_analytic_rust(int Q, const double mu, const double lambda,
                     double *dXdx_init, double *dudX, double **stored_values,
                     double *f1);
void df_analytic_rust(int Q, const double mu, const double lambda,
                      double *ddudX, double **stored_values, double *df);

// Functions defined in rust
void compute_f_analytic(int Q, const double mu, const double lambda,
                        double *dXdx_init, double *dudX, int num_comp_stored,
                        double *stored_values, double *f1);
void compute_df_analytic(int Q, const double mu, const double lambda,
                         double *ddudX, int num_comp_stored,
                         double *stored_values, double *df);

#ifdef __cplusplus
}
#endif

#endif // ANALYTIC_RUST_H
