#ifndef TAPENADE_H
#define TAPENADE_H

#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NUM_COMPONENTS_STORED_TAPENADE 15

typedef struct {
    double lambda;
    double mu;
    double *stored;
} TapenadeContext;

void init_tapenade(void *ctx);
void free_tapenade(void *ctx);
void f_tapenade(void *ctx, const double dXdx_initial[3][3], const double dudX[3][3], double f1[3][3]);
void df_tapenade(void *ctx, const double ddudX[3][3], double df1[3][3]);

double  StrainEnergy_NeoHookeanCurrentAD_Tapenade(double e_sym[6], double lambda, double mu);
void Kirchhofftau_sym_NeoHookean_AD_Tapenade(const double lambda, const double mu, double e_sym[6], double tau_sym[6]);
void dtau_fwd_Tapenade(const double lambda, const double mu, double e_sym[6], double de_sym[6], double tau_sym[6], double dtau_sym[6]);

#ifdef __cplusplus
}
#endif

#endif // TAPENADE_H
