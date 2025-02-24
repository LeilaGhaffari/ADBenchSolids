#ifndef ANALYTIC_C_H
#define ANALYTIC_C_H

#include "../utils.h"

#define NUM_COMPONENTS_STORED_ANALYTIC_C 16

void init_data_analytic_c(double **stored_values, int Q, int *num_comp);
void free_data_analytic_c(double **stored_values);
void f_analytic_c(int Q, const double mu, const double lambda, double *dXdx_init, double *dudX, double **stored_values, double *f1);
void df_analytic_c(int Q, const double mu, const double lambda, double *ddudX, double **stored_values, double *df);

#endif // ANALYTIC_C_H
