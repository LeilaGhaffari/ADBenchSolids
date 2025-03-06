#include "../../include/ad-tools/enzyme-c.h"

void init_data_enzyme_c(double **stored_values, int Q, int *num_comp) {
  *stored_values =
      (double *)malloc(Q * NUM_COMPONENTS_STORED_ENZYME_C * sizeof(double));
  *num_comp = NUM_COMPONENTS_STORED_ENZYME_C;
}

void free_data_enzyme_c(double **stored_values) {
  if (*stored_values != NULL) {
    free(*stored_values);
    *stored_values = NULL;
  }
}

double StrainEnergy_NeoHookeanCurrentAD_Enzyme(double e_sym[6], double lambda,
                                               double mu) {
  double e2_sym[6];

  // J and log(J)
  for (int i = 0; i < 6; i++)
    e2_sym[i] = 2 * e_sym[i];
  const double detbm1 = MatDetAM1Symmetric(e2_sym);
  const double J = sqrt(detbm1 + 1);
  const double logJ = Log1pSeries(detbm1) / 2.;

  // trace(e)
  const double trace_e = MatTraceSymmetric(e_sym);

  return lambda * (J * J - 1) / 4 - lambda * logJ / 2 + mu * (-logJ + trace_e);
}

void *__enzyme_function_like[2] = {(void *)Log1pSeries, (void *)"log1p"};

void Kirchhofftau_sym_NeoHookean_AD_Enzyme(const double lambda, const double mu,
                                           double e_sym[6], double tau_sym[6]) {
  double dPsi_sym[6] = {0.}, b_sym[6], dPsi[3][3], b[3][3], tau[3][3];

  // dPsi / de
  __enzyme_autodiff((void *)StrainEnergy_NeoHookeanCurrentAD_Enzyme, e_sym,
                    dPsi_sym, enzyme_const, lambda, enzyme_const, mu);
  for (int i = 3; i < 6; i++)
    dPsi_sym[i] /= 2.;

  // b = 2 e + I
  for (int j = 0; j < 6; j++)
    b_sym[j] = 2 * e_sym[j] + (j < 3);

  // tau = (dPsi / de) b
  SymmetricMatUnpack(dPsi_sym, dPsi);
  SymmetricMatUnpack(b_sym, b);
  MatMatMult(1., dPsi, b, tau);
  SymmetricMatPack(tau, tau_sym);
}

void dtau_fwd_Enzyme(const double lambda, const double mu, double e_sym[6],
                     double de_sym[6], double tau_sym[6], double dtau_sym[6]) {
  __enzyme_fwddiff((void *)Kirchhofftau_sym_NeoHookean_AD_Enzyme, enzyme_const,
                   lambda, enzyme_const, mu, e_sym, de_sym, tau_sym, dtau_sym);
}

// Residual Evaluation
void f_enzyme_c(int Q, const double mu, const double lambda, double *dXdx_init,
                double *dudX, double **stored_values, double *f1) {
  BenchPragmaSIMD for (int i = 0; i < Q; i++) {
    double Grad_u[3][3], F_inv[3][3], tau_sym[6], tau[3][3], dXdx[3][3],
        e_sym[6];
    // Pack input data
    double dudX_loc[3][3], dXdx_init_loc[3][3];
    QDataPackMat(i, Q, dXdx_init, dXdx_init_loc);
    QDataPackMat(i, Q, dudX, dudX_loc);

    MatMatMult(1.0, dudX_loc, dXdx_init_loc, Grad_u);
    double F[3][3];
    DeformationGradient(Grad_u, F);
    const double Jm1 = MatDetAM1(Grad_u);
    const double detF = Jm1 + 1.0;
    MatInverse(F, detF, F_inv);
    MatMatMult(1.0, dXdx_init_loc, F_inv, dXdx);
    GreenEulerStrain(Grad_u, e_sym);
    Kirchhofftau_sym_NeoHookean_AD_Enzyme(lambda, mu, e_sym, tau_sym);
    SymmetricMatUnpack(tau_sym, tau);
    QDataUnpackMat(i, Q, tau, f1);
    // Store
    StoredValuesPack(Q, i, 0, 9, (double *)dXdx, stored_values);
    StoredValuesPack(Q, i, 9, 6, (double *)e_sym, stored_values);
  }
}

// Jacobian Evaluation
void df_enzyme_c(int Q, const double mu, const double lambda, double *ddudX,
                 double **stored_values, double *df) {
  BenchPragmaSIMD for (int i = 0; i < Q; i++) {
    double grad_du[3][3], b_sym[6], b[3][3], de_sym[6], tau_sym[6], dtau_sym[6],
        tau[3][3], dtau[3][3], tau_grad_du[3][3], dXdx[3][3], e_sym[6],
        df_mat[3][3];
    // Unpack stored values
    StoredValuesUnpack(Q, i, 0, 9, stored_values, (double *)dXdx);
    StoredValuesUnpack(Q, i, 9, 6, stored_values, (double *)e_sym);
    // Pack input data
    double ddudX_loc[3][3];
    QDataPackMat(i, Q, ddudX, ddudX_loc);
    MatMatMult(1.0, ddudX_loc, dXdx, grad_du);
    for (int j = 0; j < 6; j++)
      b_sym[j] = 2 * e_sym[j] + (j < 3);
    SymmetricMatUnpack(b_sym, b);
    GreenEulerStrain_fwd(grad_du, b, de_sym);
    dtau_fwd_Enzyme(lambda, mu, e_sym, de_sym, tau_sym, dtau_sym);
    SymmetricMatUnpack(tau_sym, tau);
    SymmetricMatUnpack(dtau_sym, dtau);
    MatMatTransposeMult(tau, grad_du, tau_grad_du);
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        df_mat[j][k] = dtau[j][k] - tau_grad_du[j][k];
      }
    }
    QDataUnpackMat(i, Q, df_mat, df);
  }
}
