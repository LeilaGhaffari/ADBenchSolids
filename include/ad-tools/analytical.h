#ifndef ANALYTIC_H
#define ANALYTIC_H

#include "../utils.h"

#define NUM_COMPONENTS_STORED_ANALYTICAL 16

void init_data_analytic(double **stored_values, int Q);
void free_data_analytic(double **stored_values);
void f_analytic(int Q, const double mu, const double lambda, double *dXdx_init, double *dudX, double **stored_values, double *f1);
void df_analytic(int Q, const double mu, const double lambda, double *ddudX, double **stored_values, double *df);

#endif // ANALYTIC_H
