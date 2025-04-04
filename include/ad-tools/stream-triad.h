#ifndef STREAM_TRIAD_H
#define STREAM_TRIAD_H

#include "../utils.h"

#define NUM_COMPONENTS_STORED_STREAM_TRIAD 0

void init_data_stream_triad(double **stored_values, int Q, int *num_comp);
void free_data_stream_triad(double **stored_values);
void f_stream_triad(int Q, const double mu, const double lambda,
                    double *dXdx_init, double *dudX, double **stored_values,
                    double *f1);
void df_stream_triad(int Q, const double mu, const double lambda, double *ddudX,
                     double **stored_values, double *df);

#endif // STREAM_TRIAD_H
