#include "../../include/ad-tools/enzyme-rust.h"

void init_data_enzyme_rust(double **stored_values, int Q, int *num_comp) {
  *stored_values =
      (double *)malloc(Q * NUM_COMPONENTS_STORED_ENZYME_RUST * sizeof(double));
  *num_comp = NUM_COMPONENTS_STORED_ENZYME_RUST;
}

void free_data_enzyme_rust(double **stored_values) {
  if (*stored_values != NULL) {
    free(*stored_values);
    *stored_values = NULL;
  }
}

// Residual Evaluation
void f_enzyme_rust(int Q, const double mu, const double lambda,
                   double *dXdx_init, double *dudX, double **stored_values,
                   double *f1) {
  compute_f_enzyme_rust(Q, mu, lambda, dXdx_init, dudX,
                        NUM_COMPONENTS_STORED_ENZYME_RUST, *stored_values, f1);
}

// Jacobian Evaluation
void df_enzyme_rust(int Q, const double mu, const double lambda, double *ddudX,
                    double **stored_values, double *df) {
  compute_df_enzyme_rust(Q, mu, lambda, ddudX,
                         NUM_COMPONENTS_STORED_ENZYME_RUST, *stored_values, df);
}
