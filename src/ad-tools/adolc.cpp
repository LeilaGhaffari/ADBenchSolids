#include "../../include/ad-tools/adolc.h"

void init_data_adolc(double **stored_values, int Q, int *num_comp) {
  *stored_values =
      (double *)malloc(Q * NUM_COMPONENTS_STORED_ADOLC * sizeof(double));
  *num_comp = NUM_COMPONENTS_STORED_ADOLC;
}

void free_data_adolc(double **stored_values) {
  if (*stored_values != NULL) {
    free(*stored_values);
    *stored_values = NULL;
  }
}

BENCH_QFUNCTION_HELPER adouble MatDetAM1Symmetric(adouble A_sym[6]) {
  return A_sym[0] * (A_sym[1] * A_sym[2] - A_sym[3] * A_sym[3]) +
         A_sym[5] * (A_sym[3] * A_sym[4] - A_sym[5] * A_sym[2]) +
         A_sym[4] * (A_sym[5] * A_sym[3] - A_sym[4] * A_sym[1]) + A_sym[0] +
         A_sym[1] + A_sym[2] + A_sym[0] * A_sym[1] + A_sym[0] * A_sym[2] +
         A_sym[1] * A_sym[2] - A_sym[5] * A_sym[5] - A_sym[4] * A_sym[4] -
         A_sym[3] * A_sym[3];
}

BENCH_QFUNCTION_HELPER adouble MatTraceSymmetric(adouble A_sym[6]) {
  return A_sym[0] + A_sym[1] + A_sym[2];
}

BENCH_QFUNCTION_HELPER adouble Log1p_stable(adouble x) {
  adouble y, z, r;

  y = 1 + x;
  z = y - 1;
  // cancels errors with IEEE arithmetic
  r = log(y) - (z - x) / y;
  return r;
}

BENCH_QFUNCTION_HELPER adouble StrainEnergy_NeoHookeanCurrentAD_ADOLC(
    adouble e_sym[6], const double lambda, const double mu) {
  adouble E2_sym[6];
  for (int i = 0; i < 6; i++)
    E2_sym[i] = 2 * e_sym[i];
  adouble detCm1 = MatDetAM1Symmetric(E2_sym);
  adouble J = sqrt(detCm1 + 1);
  adouble logJ = Log1p_stable(detCm1) / 2.;
  adouble traceE = MatTraceSymmetric(e_sym);
  return lambda * (J * J - 1) / 4 - lambda * logJ / 2 + mu * (-logJ + traceE);
}

BENCH_QFUNCTION_HELPER void ComputeGradPsi(double grad[6], double e_sym[6],
                                           const double lambda,
                                           const double mu) {
  double Fp[1];
  adouble ea[6], Fa[1];
  int tag = 1;
  trace_on(tag);
  for (int i = 0; i < 6; i++)
    ea[i] <<= e_sym[i];
  Fa[0] = StrainEnergy_NeoHookeanCurrentAD_ADOLC(ea, lambda, mu);
  Fa[0] >>= Fp[0];
  trace_off();
  gradient(tag, 6, e_sym, grad);
  for (int i = 0; i < 6; i++)
    if (i > 2)
      grad[i] /= 2.;
}

BENCH_QFUNCTION_HELPER void ComputeHessianPsi(double hess[6][6],
                                              double e_sym[6],
                                              const double lambda,
                                              const double mu) {
  adouble ea[6], Fa[1];
  double Fp[1], buf[21];
  double *H[6] = {&buf[0], &buf[1], &buf[3], &buf[6], &buf[10], &buf[15]};
  int tag = 1;
  trace_on(tag);
  for (int i = 0; i < 6; i++)
    ea[i] <<= e_sym[i];
  Fa[0] = StrainEnergy_NeoHookeanCurrentAD_ADOLC(ea, lambda, mu);
  Fa[0] >>= Fp[0];
  trace_off();
  hessian(tag, 6, e_sym, H);

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < i + 1; j++) {
      hess[i][j] = H[i][j];
      if (i != j)
        hess[j][i] = hess[i][j];
    }
  }
  for (int i = 3; i < 6; i++)
    for (int j = 0; j < 6; j++)
      hess[i][j] /= 2.;
}

void f_adolc(int Q, const double mu, const double lambda, double *dXdx_init,
             double *dudX, double **stored_values, double *f1) {
  BenchPragmaSIMD for (int i = 0; i < Q; i++) {
    double Grad_u[3][3], F_inv[3][3], tau_sym[6], dXdx[3][3], e_sym[6];
    // Pack input data
    double dudX_loc[3][3], dXdx_init_loc[3][3];
    QDataPackMat(i, Q, dXdx_init, dXdx_init_loc);
    QDataPackMat(i, Q, dudX, dudX_loc);

    // Grad_u = du/dx_initial = du/dX * dX/dx_initial
    MatMatMult(1.0, dudX_loc, dXdx_init_loc, Grad_u);
    // Compute the Deformation Gradient : F = I + Grad_u
    double F[3][3];
    DeformationGradient(Grad_u, F);
    // Compute F^{-1}
    const double Jm1 = MatDetAM1(Grad_u);
    const double detF = Jm1 + 1.;
    MatInverse(F, detF, F_inv);
    // dXdx = dX/dx = dX/dx_initial * F^{-1}
    MatMatMult(1.0, dXdx_init_loc, F_inv, dXdx);
    // Compute Green Euler Strain tensor (e)
    GreenEulerStrain(Grad_u, e_sym);
    // Compute gradient of Psi (dPsi/de) with ADOL-C
    double gradPsi_sym[6] = {0.};
    ComputeGradPsi(gradPsi_sym, e_sym, lambda, mu);
    // b = 2 e + I
    double b_sym[6], b[3][3];
    for (int j = 0; j < 6; j++)
      b_sym[j] = 2. * e_sym[j] + (j < 3);
    SymmetricMatUnpack(b_sym, b);
    // tau = dPsi/de * b
    double tau[3][3], gradPsi[3][3];
    SymmetricMatUnpack(gradPsi_sym, gradPsi);
    MatMatMult(1.0, gradPsi, b, tau);
    SymmetricMatPack(tau, tau_sym);
    SymmetricMatUnpack(tau_sym, tau);
    QDataUnpackMat(i, Q, tau, f1);
    // Store
    StoredValuesPack(Q, i, 0, 9, (double *)dXdx, stored_values);
    StoredValuesPack(Q, i, 9, 6, (double *)e_sym, stored_values);
    StoredValuesPack(Q, i, 15, 6, (double *)gradPsi_sym, stored_values);
  }
}

void df_adolc(int Q, const double mu, const double lambda, double *ddudX,
              double **stored_values, double *df) {
  BenchPragmaSIMD for (int i = 0; i < Q; i++) {
    double grad_du[3][3], b_sym[6], b[3][3], de_sym[6], dtau_sym[6], tau[3][3],
        dtau[3][3], tau_grad_du[3][3], dXdx[3][3], e_sym[6], gradPsi_sym[6],
        df_mat[3][3];
    // Unpack stored values
    StoredValuesUnpack(Q, i, 0, 9, stored_values, (double *)dXdx);
    StoredValuesUnpack(Q, i, 9, 6, stored_values, (double *)e_sym);
    StoredValuesUnpack(Q, i, 15, 6, stored_values, (double *)gradPsi_sym);
    // Pack input data
    double ddudX_loc[3][3];
    QDataPackMat(i, Q, ddudX, ddudX_loc);
    // Compute grad_du = ddu/dX * dX/dx
    // X is ref coordinate [-1,1]^3; x is physical coordinate in current
    // configuration
    MatMatMult(1.0, ddudX_loc, dXdx, grad_du);
    // Compute b = 2 e + I
    for (int j = 0; j < 6; j++)
      b_sym[j] = 2 * e_sym[j] + (j < 3);
    SymmetricMatUnpack(b_sym, b);
    // Compute de = db / 2 = (grad_du b + b (grad_du)^T) / 2
    double de[3][3];
    GreenEulerStrain_fwd(grad_du, b, de_sym);
    SymmetricMatUnpack(de_sym, de);
    // ADOL-C
    double hessPsi_curr[6][6] = {{0.}};
    ComputeHessianPsi(hessPsi_curr, e_sym, lambda, mu); // hessian: d2Psi/de2
    //---------------------------------------------------
    // dtau = (hessPsi : de) b + 2 gradPsi (I_4 : de)
    //      = dGradPsi b + 2 gradPsi de
    //      = dtau_1 + dtau_2
    //---------------------------------------------------
    // dGradPsi = hessPsi : de
    double dGradPsi[3][3], dGradPsi_sym[6] = {0.};
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 6; j++) {
        dGradPsi_sym[i] += hessPsi_curr[i][j] * de_sym[j];
      }
    }
    SymmetricMatUnpack(dGradPsi_sym, dGradPsi);
    // dtau_1 = dGradPsi b
    double dtau_1[3][3];
    MatMatMult(1.0, dGradPsi, b, dtau_1);
    // dtau_2 = 2 gradPsi de
    double dtau_2[3][3], gradPsi[3][3];
    SymmetricMatUnpack(gradPsi_sym, gradPsi);
    MatMatMult(2., gradPsi, de, dtau_2);
    // dtau = dtau_1 + dtau_2
    MatMatAdd(1., dtau_1, 1., dtau_2, dtau);
    SymmetricMatPack(dtau, dtau_sym);
    //---------------------------------------------------
    // tau = dPsi/de * b
    MatMatMult(1.0, gradPsi, b, tau);
    // Compute tau_grad_du = tau * grad_du^T
    MatMatTransposeMult(tau, grad_du, tau_grad_du);
    // Compute df = dtau - tau * grad_du^T
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        df_mat[j][k] = dtau[j][k] - tau_grad_du[j][k];
      }
    }
    QDataUnpackMat(i, Q, df_mat, df);
  }
}
