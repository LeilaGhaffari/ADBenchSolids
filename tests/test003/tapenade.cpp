#include "tapenade.h"
#include <math.h>
#include <stdlib.h>

// Residual Evaluation
void f_tapenade(void *ctx, const double dXdx_initial[3][3],
                const double dudX[3][3], double f1[3][3]) {
  // Context
  TapenadeContext *context = static_cast<TapenadeContext *>(ctx);
  const double mu = context->mu;
  const double lambda = context->lambda;
  double *stored_values = context->stored;

  double Grad_u[3][3], F_inv[3][3], tau_sym[6], dXdx[3][3], e_sym[6];

  MatMatMult(1.0, dudX, dXdx_initial, Grad_u);

  double F[3][3];
  DeformationGradient(Grad_u, F);

  const double Jm1 = MatDetAM1(Grad_u);
  const double detF = Jm1 + 1.0;
  MatInverse(F, detF, F_inv);

  MatMatMult(1.0, dXdx_initial, F_inv, dXdx);

  GreenEulerStrain(Grad_u, e_sym);

  kirchhofftau_tapenade_(lambda, mu, e_sym, tau_sym);
  SymmetricMatUnpack(tau_sym, f1);

  // ------------------------------------------------------------------------
  // Store
  // ------------------------------------------------------------------------
  StoredValuesPack(0, 9, (double *)dXdx, stored_values);
  StoredValuesPack(9, 6, (double *)e_sym, stored_values);

  // ------------------------------------------------------------------------
  // More info
  // ------------------------------------------------------------------------
  printf("\n\ne =");
  for (int i = 0; i < 6; i++)
    printf("\n\t%.12lf", e_sym[i]);
  printf("\n\nStrain Energy from e = ");
  printf(" %.12lf", StrainEnergy(e_sym, context->lambda, context->mu));
  printf("\n\ntau =");
  for (int i = 0; i < 6; i++)
    printf("\n\t%.12lf", tau_sym[i]);
}

// Jacobian Evaluation
void df_tapenade(void *ctx, const double ddudX[3][3], double df1[3][3]) {
  // Context
  TapenadeContext *context = static_cast<TapenadeContext *>(ctx);
  const double mu = context->mu;
  const double lambda = context->lambda;
  const double *stored_values = context->stored;

  double grad_du[3][3], b_sym[6], b[3][3], de_sym[6], tau_sym[6], dtau_sym[6],
      tau[3][3], dtau[3][3], tau_grad_du[3][3], dXdx[3][3], e_sym[6];
  // ------------------------------------------------------------------------
  // Unpack stored values
  // ------------------------------------------------------------------------
  StoredValuesUnpack(0, 9, stored_values, (double *)dXdx);
  StoredValuesUnpack(9, 6, stored_values, (double *)e_sym);

  MatMatMult(1.0, ddudX, dXdx, grad_du);

  for (int j = 0; j < 6; j++)
    b_sym[j] = 2 * e_sym[j] + (j < 3);
  SymmetricMatUnpack(b_sym, b);

  GreenEulerStrain_fwd(grad_du, b, de_sym);

  dtau_fwd_tapenade_(lambda, mu, e_sym, de_sym, tau_sym, dtau_sym);

  SymmetricMatUnpack(tau_sym, tau);
  SymmetricMatUnpack(dtau_sym, dtau);

  MatMatTransposeMult(tau, grad_du, tau_grad_du);
  for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {
      df1[j][k] = dtau[j][k] - tau_grad_du[j][k];
    }
  }

  // ------------------------------------------------------------------------
  // More info
  // ------------------------------------------------------------------------
  // printf("\n\nde =");
  // for (int i=0; i<6; i++) printf("\n\t%.12lf", de_sym[i]);
  // printf("\n\ndtau =");
  // for (int i=0; i<6; i++) printf("\n\t%.12lf", dtau_sym[i]);
  printf("\n\n");
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++)
      printf("\t%.12lf", df1[i][j]);
    printf("\n");
  }
  printf("\n\n");
}
