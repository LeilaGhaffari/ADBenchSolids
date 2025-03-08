#ifndef ENZYME_C_H
#define ENZYME_C_H

#include "../utils.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NUM_COMPONENTS_STORED_ENZYME_C 15

void init_data_enzyme_c(double **stored_values, int Q, int *num_comp);
void free_data_enzyme_c(double **stored_values);
void f_enzyme_c(int Q, const double mu, const double lambda, double *dXdx_init,
                double *dudX, double **stored_values, double *f1);
void df_enzyme_c(int Q, const double mu, const double lambda, double *ddudX,
                 double **stored_values, double *df1);

#ifdef __cplusplus
}
#endif

#endif // ENZYME_C_H
