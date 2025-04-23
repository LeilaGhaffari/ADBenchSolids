#ifndef ADOLC_H
#define ADOLC_H

#include "../utils.h"
#include <adolc/adolc.h>

#define NUM_COMPONENTS_STORED_ADOLC 15

void init_data_adolc(double **stored_values, int Q, int *num_comp);
void free_data_adolc(double **stored_values);
void f_adolc(int Q, const double mu, const double lambda, double *dXdx_init,
             double *dudX, double **stored_values, double *f1);
void df_adolc(int Q, const double mu, const double lambda, double *ddudX,
              double **stored_values, double *df);

#endif // ADOLC_H
