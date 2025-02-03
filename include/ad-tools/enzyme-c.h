#ifndef ENZYME_C_H
#define ENZYME_C_H

#include "../utils.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NUM_COMPONENTS_STORED_ENZYME_C 15

void init_data_enzyme_c(double **stored_values, int Q);
void free_data_enzyme_c(double **stored_values);
void f_enzyme_c(int Q, const double mu, const double lambda, double *dXdx_init, double *dudX, double **stored_values, double *f1);
void df_enzyme_c(int Q, const double mu, const double lambda, double *ddudX, double **stored_values, double *df1);

void __enzyme_autodiff(void *, ...);
void __enzyme_fwddiff(void *, ...);
extern int enzyme_const;
double  StrainEnergy_NeoHookeanCurrentAD_Enzyme(double e_sym[6], double lambda, double mu);
void Kirchhofftau_sym_NeoHookean_AD_Enzyme(const double lambda, const double mu, double e_sym[6], double tau_sym[6]);
void dtau_fwd_Enzyme(const double lambda, const double mu, double e_sym[6], double de_sym[6], double tau_sym[6], double dtau_sym[6]);

#ifdef __cplusplus
}
#endif

#endif // ENZYME_C_H
