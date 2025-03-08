#include "../../include/ad-tools/tapenade.h"

void init_data_tapenade(double **stored_values, int Q, int *num_comp) {
  *stored_values =
      (double *)malloc(Q * NUM_COMPONENTS_STORED_TAPENADE * sizeof(double));
  *num_comp = NUM_COMPONENTS_STORED_TAPENADE;
}

void free_data_tapenade(double **stored_values) {
  if (*stored_values != NULL) {
    free(*stored_values);
    *stored_values = NULL;
  }
}

BENCH_QFUNCTION_HELPER void tau_sym_ad(const double e_sym[6],
                                       const double lambda, const double mu,
                                       double tau_sym[6]) {
  double lambdab = 0., mub = 0., energy, energyb = 1., grad_psi_sym[6] = {0.};
  for (int i = 0; i < 6; i++)
    grad_psi_sym[i] = 0.;
  compute_grad_psi_tapenade(e_sym, grad_psi_sym, lambda, &lambdab, mu, &mub,
                            &energy, &energyb);
  for (int i = 3; i < 6; i++)
    grad_psi_sym[i] /= 2.;

  // b = 2 e + I
  double b_sym[6];
  for (int j = 0; j < 6; j++)
    b_sym[j] = 2 * e_sym[j] + (j < 3);

  // tau = (dPsi / de) b
  double grad_psi[3][3], b[3][3], tau[3][3];
  SymmetricMatUnpack_t(grad_psi_sym, grad_psi);
  SymmetricMatUnpack_t(b_sym, b);
  MatMatMult_t(1., grad_psi, b, tau);
  SymmetricMatPack_t(tau, tau_sym);
}

BENCH_QFUNCTION_HELPER void
dtau_sym_fwd(const double e_sym[6], const double de_sym[6], const double lambda,
             const double mu, double tau_sym[6], double dtau_sym[6]) {
  const double lambdad = 0., mud = 0;
  compute_dtau_sym_fwd_tapenade(e_sym, de_sym, lambda, lambdad, mu, mud,
                                tau_sym, dtau_sym);
}

// Residual Evaluation
void f_tapenade(int Q, const double mu, const double lambda, double *dXdx_init,
                double *dudX, double **stored_values, double *f1) {
  BenchPragmaSIMD for (int i = 0; i < Q; i++) {
    double Grad_u[3][3], F_inv[3][3], tau_sym[6], tau[3][3], dXdx[3][3],
        e_sym[6], F[3][3], dudX_loc[3][3], dXdx_init_loc[3][3];
    // Pack input data
    QDataPackMat(i, Q, dXdx_init, dXdx_init_loc);
    QDataPackMat(i, Q, dudX, dudX_loc);

    MatMatMult(1.0, dudX_loc, dXdx_init_loc, Grad_u);
    DeformationGradient(Grad_u, F);
    const double Jm1 = MatDetAM1(Grad_u);
    const double detF = Jm1 + 1.0;
    MatInverse(F, detF, F_inv);
    MatMatMult(1.0, dXdx_init_loc, F_inv, dXdx);
    GreenEulerStrain(Grad_u, e_sym);
    tau_sym_ad(e_sym, lambda, mu, tau_sym);
    SymmetricMatUnpack(tau_sym, tau);
    QDataUnpackMat(i, Q, tau, f1);

    // Store
    StoredValuesPack(Q, i, 0, 9, (double *)dXdx, stored_values);
    StoredValuesPack(Q, i, 9, 6, (double *)e_sym, stored_values);
  }
}

// Jacobian Evaluation
void df_tapenade(int Q, const double mu, const double lambda, double *ddudX,
                 double **stored_values, double *df) {
  BenchPragmaSIMD for (int i = 0; i < Q; i++) {
    double grad_du[3][3], b_sym[6], b[3][3], de_sym[6], tau_sym[6], dtau_sym[6],
        tau[3][3], dtau[3][3], tau_grad_du[3][3], dXdx[3][3], e_sym[6],
        df_mat[3][3], ddudX_loc[3][3];
    // Unpack stored values
    StoredValuesUnpack(Q, i, 0, 9, stored_values, (double *)dXdx);
    StoredValuesUnpack(Q, i, 9, 6, stored_values, (double *)e_sym);

    // Pack input data
    QDataPackMat(i, Q, ddudX, ddudX_loc);

    MatMatMult(1.0, ddudX_loc, dXdx, grad_du);
    for (int j = 0; j < 6; j++)
      b_sym[j] = 2 * e_sym[j] + (j < 3);
    SymmetricMatUnpack(b_sym, b);
    GreenEulerStrain_fwd(grad_du, b, de_sym);
    dtau_sym_fwd(e_sym, de_sym, lambda, mu, tau_sym, dtau_sym);
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
