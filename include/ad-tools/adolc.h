#ifndef ADOLC_H
#define ADOLC_H

#include "../utils.h"
#include <adolc/adolc.h>

#define NUM_COMPONENTS_STORED_ADOLC 21

void init_data_adolc(double **stored_values, int Q, int *num_comp);
void free_data_adolc(double **stored_values);
void f_adolc(int Q, const double mu, const double lambda, double *dXdx_init,
             double *dudX, double **stored_values, double *f1);
void df_adolc(int Q, const double mu, const double lambda, double *ddudX,
              double **stored_values, double *df);

adouble MatDetAM1Symmetric(adouble A_sym[6]);
adouble MatTraceSymmetric(adouble A_sym[6]);
adouble Log1pSeries(adouble x);
adouble StrainEnergy_NeoHookeanCurrentAD_ADOLC(adouble e_sym[6],
                                               const double lambda,
                                               const double mu);
void ComputeGradPsi(double grad[6], double e_sym[6], const double lambda,
                    const double mu);
void ComputeHessianPsi(double hess[6][6], double e_sym[6], const double lambda,
                       const double mu);

#endif // ADOLC_H
