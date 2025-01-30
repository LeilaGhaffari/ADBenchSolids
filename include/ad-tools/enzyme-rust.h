#ifndef ENZYME_RUST_H
#define ENZYME_RUST_H

#include "../utils.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NUM_COMPONENTS_STORED_ENZYME_RUST 15

void init_data_enzyme_rust(double **stored_values, int Q);
void free_data_enzyme_rust(double **stored_values);
void f_enzyme_rust(int Q, const double mu, const double lambda, double *dXdx_init, double *dudX, double **stored_values, double *f1);
void df_enzyme_rust(int Q, const double mu, const double lambda, double *ddudX, double **stored_values, double *df);

void compute_stress(double lambda, double mu, const double e_sym[6], double tau_sym[6]);

#ifdef __cplusplus
}
#endif

#endif // ENZYME_RUST_H
