#ifndef STREAM_RESIDUAL_H
#define STREAM_RESIDUAL_H

#include "../utils.h"

#define NUM_COMPONENTS_STORED_STREAM_RESIDUAL 15

void init_data_stream_residual(double **stored_values, int Q, int *num_comp);
void free_data_stream_residual(double **stored_values);
void f_stream_residual(int Q, const double mu, const double lambda,
                       double *dXdx_init, double *dudX, double **stored_values,
                       double *f1);
void df_stream_residual(int Q, const double mu, const double lambda,
                        double *ddudX, double **stored_values, double *df);

#endif // STREAM_RESIDUAL_H
