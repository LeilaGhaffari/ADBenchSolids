#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline int QDataPackMat(int i, double *in, double out[3][3]) {
  int idx_start = 9*i;
  out[0][0] = in[idx_start + 0];
  out[0][1] = in[idx_start + 1];
  out[0][2] = in[idx_start + 2];
  out[1][0] = in[idx_start + 3];
  out[1][1] = in[idx_start + 4];
  out[1][2] = in[idx_start + 5];
  out[2][0] = in[idx_start + 6];
  out[2][1] = in[idx_start + 7];
  out[2][2] = in[idx_start + 8];
  return 0;
}

static inline int QDataUnpackMat(int i, double in[3][3], double *out) {
  int idx_start = 9*i;
  out[idx_start + 0] = in[0][0];
  out[idx_start + 1] = in[0][1];
  out[idx_start + 2] = in[0][2];
  out[idx_start + 3] = in[1][0];
  out[idx_start + 4] = in[1][1];
  out[idx_start + 5] = in[1][2];
  out[idx_start + 6] = in[2][0];
  out[idx_start + 7] = in[2][1];
  out[idx_start + 8] = in[2][2];
  return 0;
}

static inline int SymmetricMatUnpack(const double sym[6], double full[3][3]) {
  full[0][0] = sym[0];
  full[0][1] = sym[5];
  full[0][2] = sym[4];
  full[1][0] = sym[5];
  full[1][1] = sym[1];
  full[1][2] = sym[3];
  full[2][0] = sym[4];
  full[2][1] = sym[3];
  full[2][2] = sym[2];
  return 0;
};

static inline int SymmetricMatPack(const double full[3][3], double sym[6]) {
  sym[0] = full[0][0];
  sym[1] = full[1][1];
  sym[2] = full[2][2];
  sym[3] = full[1][2];
  sym[4] = full[0][2];
  sym[5] = full[0][1];
  return 0;
};

static inline int MatInverse(const double A[3][3], double det_A, double A_inv[3][3]) {
  double B[3][3] = {
      {A[1][1] * A[2][2] - A[1][2] * A[2][1], A[0][2] * A[2][1] - A[0][1] * A[2][2], A[0][1] * A[1][2] - A[0][2] * A[1][1]},
      {A[1][2] * A[2][0] - A[1][0] * A[2][2], A[0][0] * A[2][2] - A[0][2] * A[2][0], A[0][2] * A[1][0] - A[0][0] * A[1][2]},
      {A[1][0] * A[2][1] - A[1][1] * A[2][0], A[0][1] * A[2][0] - A[0][0] * A[2][1], A[0][0] * A[1][1] - A[0][1] * A[1][0]},
  };
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A_inv[i][j] = B[i][j] / (det_A);
    }
  }
  return 0;
};

static inline int MatMatMult(double alpha, const double A[3][3], const double B[3][3], double C[3][3]) {
  for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {
      C[j][k] = 0;
      for (int m = 0; m < 3; m++) {
        C[j][k] += alpha * A[j][m] * B[m][k];
      }
    }
  }
  return 0;
};

static inline int MatComputeInverseSymmetric(const double A[3][3], const double det_A, double A_inv[6]) {
  // Compute A^(-1) : A-Inverse
  double B[6] = {
      A[1][1] * A[2][2] - A[1][2] * A[2][1],
      A[0][0] * A[2][2] - A[0][2] * A[2][0],
      A[0][0] * A[1][1] - A[0][1] * A[1][0],
      A[0][2] * A[1][0] - A[0][0] * A[1][2],
      A[0][1] * A[1][2] - A[0][2] * A[1][1],
      A[0][2] * A[2][1] - A[0][1] * A[2][2]
  };
  for (int m = 0; m < 6; m++) {
    A_inv[m] = B[m] / (det_A);
  }
  return 0;
};

static inline double MatDetAM1Symmetric(const double A_sym[6]) {
  return A_sym[0] * (A_sym[1] * A_sym[2] - A_sym[3] * A_sym[3]) +
         A_sym[5] * (A_sym[3] * A_sym[4] - A_sym[5] * A_sym[2]) +
         A_sym[4] * (A_sym[5] * A_sym[3] - A_sym[4] * A_sym[1]) +
         A_sym[0] + A_sym[1] + A_sym[2] +
         A_sym[0] * A_sym[1] + A_sym[0] * A_sym[2] + A_sym[1] * A_sym[2] -
         A_sym[5] * A_sym[5] - A_sym[4] * A_sym[4] - A_sym[3] * A_sym[3];
};

static inline double MatTraceSymmetric(const double A_sym[6]) { return A_sym[0] + A_sym[1] + A_sym[2]; };

static inline double MatDetAM1(const double A[3][3]) {
    return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) +
           A[0][1] * (A[1][2] * A[2][0] - A[1][0] * A[2][2]) +
           A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]) +
           A[0][0] + A[1][1] + A[2][2] +
           A[0][0] * A[1][1] + A[0][0] * A[2][2] + A[1][1] * A[2][2] -
           A[0][1] * A[1][0] - A[0][2] * A[2][0] - A[1][2] * A[2][1];
};

static inline void DeformationGradient(double grad_u[3][3], double F[3][3]) {
  F[0][0] = grad_u[0][0] + 1.;
  F[0][1] = grad_u[0][1];
  F[0][2] = grad_u[0][2];
  F[1][0] = grad_u[1][0];
  F[1][1] = grad_u[1][1] + 1.;
  F[1][2] = grad_u[1][2];
  F[2][0] = grad_u[2][0];
  F[2][1] = grad_u[2][1];
  F[2][2] = grad_u[2][2] + 1.;
};

static inline double Log1pSeries(double x) {
    double sum = 0;
    double y = x / (2. + x);
    double y2 = y*y;
    sum += y;
    for (int i=0; i<5; i++) {
      y *= y2;
      sum += y / (2*i + 3);
    }
    return 2 * sum;
};

static inline int LinearStrain(const double grad_u[3][3], double e_sym[6]) {
  e_sym[0] = grad_u[0][0];
  e_sym[1] = grad_u[1][1];
  e_sym[2] = grad_u[2][2];
  e_sym[3] = (grad_u[1][2] + grad_u[2][1]) / 2.;
  e_sym[4] = (grad_u[0][2] + grad_u[2][0]) / 2.;
  e_sym[5] = (grad_u[0][1] + grad_u[1][0]) / 2.;
  return 0;
};

static inline int GreenEulerStrain(const double grad_u[3][3], double e_sym[6]) {
  const int ind_j[6] = {0, 1, 2, 1, 0, 0}, ind_k[6] = {0, 1, 2, 2, 2, 1};
  LinearStrain(grad_u, e_sym);
  // Add (grad_u * grad_u^T)/2 term to the linear part of e_sym
  for (int m = 0; m < 6; m++) {
    for (int n = 0; n < 3; n++) {
      e_sym[m] += (grad_u[ind_j[m]][n] * grad_u[ind_k[m]][n]) * 0.5;
    }
  }
  return 0;
};

static inline int GreenEulerStrain_fwd(const double grad_du[3][3], const double b[3][3], double de_sym[6]) {
  const int ind_j[6] = {0, 1, 2, 1, 0, 0}, ind_k[6] = {0, 1, 2, 2, 2, 1};
  for (int m = 0; m < 6; m++) {
    de_sym[m] = 0;
    for (int n = 0; n < 3; n++) {
      de_sym[m] += (grad_du[ind_j[m]][n] * b[n][ind_k[m]] + b[ind_j[m]][n] * grad_du[ind_k[m]][n]) / 2.;
    }
  }
  return 0;
};

static inline double StrainEnergy(double e_sym[6], const double lambda, const double mu) {
  double e2_sym[6];
  // J and log(J)
  for (int i = 0; i < 6; i++) e2_sym[i] = 2 * e_sym[i];
  const double detbm1 = MatDetAM1Symmetric(e2_sym);
  const double J      = sqrt(detbm1 + 1);
  const double logJ   = Log1pSeries(detbm1) / 2.;
  // trace(e)
  const double trace_e = MatTraceSymmetric(e_sym);
  return lambda * (J * J - 1) / 4 - lambda * logJ / 2 + mu * (-logJ + trace_e);
};

static inline int MatMatTransposeMult(const double A[3][3], const double B[3][3], double C[3][3]) {
  for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {
      C[j][k] = 0;
      for (int m = 0; m < 3; m++) {
        C[j][k] += A[j][m] * B[k][m];
      }
    }
  }
  return 0;
};

static inline void PullBack_symmetric(double Grad_u[3][3], double a_sym[6], double A_sym[6]) {
  // F = I + Grad_u
  const double F[3][3] = {
    {Grad_u[0][0] + 1, Grad_u[0][1],     Grad_u[0][2]    },
    {Grad_u[1][0],     Grad_u[1][1] + 1, Grad_u[1][2]    },
    {Grad_u[2][0],     Grad_u[2][1],     Grad_u[2][2] + 1}
  };
  // F^{-1}
  double F_inv[3][3];
  const double Jm1  = MatDetAM1(Grad_u);
  const double detF = Jm1 + 1.;
  MatInverse(F, detF, F_inv);
  // A = F_inv * a * F_inv_T (pull-back)
  double F_inv_a[3][3], a[3][3], A[3][3];
  SymmetricMatUnpack(a_sym, a);
  MatMatMult(1., F_inv, a, F_inv_a);
  MatMatTransposeMult(F_inv_a, F_inv, A);
  SymmetricMatPack(A, A_sym);
};

static inline void dPullBack_symmetric(double F_inv[3][3], double dF[3][3], double a_sym[6], double da_sym[6], double dA_sym[6]) {
  // F = I + Grad_u => dF = Grad_du
  // F_inv * F = I => dF_inv * F + F_inv * dF = 0 => dF_inv = -F_inv * dF * F_inv
  // A = F_inv * a * F_inv^T => da = F_inv da F_inv^T + dF_inv a F_inv^T + F_inv a dF_inv^T
  double a[3][3], da[3][3];
  SymmetricMatUnpack(da_sym, da);
  SymmetricMatUnpack(a_sym, a);

  // dF_inv = -F_inv * dF * F_inv
  double F_inv_dF[3][3], dF_inv[3][3];
  MatMatMult(-1., F_inv, dF, F_inv_dF);
  MatMatMult(1., F_inv_dF, F_inv, dF_inv);

  // F_inv da F_inv^T
  double F_inv_da[3][3], F_inv_da_F_inv_T[3][3];
  MatMatMult(1., F_inv, da, F_inv_da);
  MatMatTransposeMult(F_inv_da, F_inv, F_inv_da_F_inv_T);

  // dF_inv a F_inv^T
  double dF_inv_a[3][3], dF_inv_a_F_inv_T[3][3];
  MatMatMult(1., dF_inv, a, dF_inv_a);
  MatMatTransposeMult(dF_inv_a, F_inv, dF_inv_a_F_inv_T);

  // F_inv a dF_inv^T
  double F_inv_a[3][3], F_inv_a_dF_inv_T[3][3];
  MatMatMult(1., F_inv, a, F_inv_a);
  MatMatTransposeMult(F_inv_a, dF_inv, F_inv_a_dF_inv_T);

  // da = F_inv da F_inv^T + dF_inv a F_inv^T + F_inv a dF_inv^T
  double dA[3][3] = {{0.}};
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      dA[i][j] = F_inv_da_F_inv_T[i][j] + dF_inv_a_F_inv_T[i][j] + F_inv_a_dF_inv_T[i][j];

  SymmetricMatPack(dA, dA_sym);
};

static inline void PushForward_symmetric(double F[3][3], double A_sym[6], double a_sym[6]) {
  // a = F * A * F^T (push-forward)
  double F_A[3][3], a[3][3], A[3][3];
  SymmetricMatUnpack(A_sym, A);
  MatMatMult(1., F, A, F_A);
  MatMatTransposeMult(F_A, F, a);
  SymmetricMatPack(a, a_sym);
};

static inline void dPushForward_symmetric(const double F[3][3], double dF[3][3], double A_sym[6], double dA_sym[6], double da_sym[6]) {
  // F = I + Grad_u => dF = Grad_du
  // a = F * A * F^T => da = F dA F^T + dF A F^T + F A dF^T
  double A[3][3], dA[3][3];
  SymmetricMatUnpack(dA_sym, dA);
  SymmetricMatUnpack(A_sym, A);

  // F dA F^T
  double F_dA[3][3], F_dA_FT[3][3];
  MatMatMult(1., F, dA, F_dA);
  MatMatTransposeMult(F_dA, F, F_dA_FT);

  // dF A F^T
  double dF_A[3][3], dF_A_FT[3][3];
  MatMatMult(1., dF, A, dF_A);
  MatMatTransposeMult(dF_A, F, dF_A_FT);

  // F A dF^T
  double F_A[3][3], F_A_dFT[3][3];
  MatMatMult(1., F, A, F_A);
  MatMatTransposeMult(F_A, dF, F_A_dFT);

  // da = F dA F^T + dF A F^T + F A dF^T
  double da[3][3] = {{0.}};
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      da[i][j] = F_dA_FT[i][j] + dF_A_FT[i][j] + F_A_dFT[i][j];

  SymmetricMatPack(da, da_sym);
};

static inline void MatTransposeMatMult(double alpha, const double A[3][3], const double B[3][3], double C[3][3]) {
  for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {
      C[j][k] = 0;
      for (int m = 0; m < 3; m++) {
        C[j][k] += alpha * A[m][j] * B[m][k];
      }
    }
  }
};

// C = alpha A + beta B for 3x3 matrices
static inline int MatMatAdd(double alpha, const double A[3][3], double beta, const double B[3][3], double C[3][3]) {
  for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {
      C[j][k] = alpha * A[j][k] + beta * B[j][k];
    }
  }
  return 0;
};

static inline int StoredValuesPack(int start, int num_comp, int num_stored_comp, int i, const double *local, double **stored) {
  start += i*num_stored_comp;
  for (int j = 0; j < num_comp; j++) (*stored)[start + j] = local[j];
  return 0;
};

static inline int StoredValuesUnpack(int start, int num_comp, int num_stored_comp, int i, double **stored, double *local) {
  start += i*num_stored_comp;
  for (int j = 0; j < num_comp; j++) local[j] = (*stored)[start + j];
  return 0;
}

/*
Compute V, J * dV/dJ, J^2 * d^2V/dJ^2 for mixed/single fields hyperelasticity.
  For isochoric single field and mixed methods we consider the volumetric energy of form
  `\psi_vol = k * V`

 `A = J^2 - 1 - 2 logJ`
 `V = A / 4`
 `J dV/dJ = (J^2 - 1) / 2`
 `J^2 d^2V/dJ^2 = (J^2 + 1) / 2`
 */
static inline int VolumetricFunctionAndDerivatives(double Jm1, double *V, double *J_dVdJ, double *J2_d2VdJ2) {
  const double Jp1  = Jm1 + 2.;
  const double logJ = Log1pSeries(Jm1);
  const double A    = Jm1 * Jp1 - 2. * logJ;

  if (V) *V = A * 0.25;
  if (J_dVdJ) *J_dVdJ = Jm1 * Jp1 * 0.5;
  if (J2_d2VdJ2) {
    const double J = Jm1 + 1.;

    *J2_d2VdJ2 = (J * J + 1) * 0.5;
  }
  return 0;
};

static inline int VolumetricKirchhoffTau(double J_dVdJ, double bulk, double *tau_vol_sym) {
  // -- [bulk * J dV/dJ]
  *tau_vol_sym = bulk * J_dVdJ;
  return 0;
};

static inline int KirchhoffTau_NeoHookean(double J_dVdJ, double lambda, double two_mu, const double e_sym[6],
                                                       double tau_sym[6]) {
  double tau_vol_sym;

  VolumetricKirchhoffTau(J_dVdJ, lambda, &tau_vol_sym);

  tau_sym[0] = two_mu * e_sym[0] + tau_vol_sym;
  tau_sym[1] = two_mu * e_sym[1] + tau_vol_sym;
  tau_sym[2] = two_mu * e_sym[2] + tau_vol_sym;
  tau_sym[3] = two_mu * e_sym[3];
  tau_sym[4] = two_mu * e_sym[4];
  tau_sym[5] = two_mu * e_sym[5];
  return 0;
};

static inline int MatMatAddSymmetric(double alpha, const double A_sym[6], double beta, const double B_sym[6],
                                                  double C_sym[6]) {
  for (int j = 0; j < 6; j++) {
    C_sym[j] = alpha * A_sym[j] + beta * B_sym[j];
  }
  return 0;
};

static inline int ScalarMatMultSymmetric(double alpha, const double A_sym[6], double B_sym[6]) {
  for (int i = 0; i < 6; i++) B_sym[i] = alpha * A_sym[i];
  return 0;
};

static inline int ComputeFdSFTransposeVolumetric(double J_dVdJ, double J2_d2VdJ2, double bulk, double trace_depsilon,
                                                              const double depsilon_sym[6], double FdSFTranspose_vol_sym[6]) {
  const double coeff_1 = (bulk * J2_d2VdJ2 + bulk * J_dVdJ) * trace_depsilon;
  const double coeff_2 = bulk * J_dVdJ;

  // F*dS_vol*F^T = -2 * coeff_2 * depsilon_sym
  ScalarMatMultSymmetric(-2. * coeff_2, depsilon_sym, FdSFTranspose_vol_sym);

  // Add coeff_1 * I to FdSFTranspose_vol_sym
  FdSFTranspose_vol_sym[0] += coeff_1;
  FdSFTranspose_vol_sym[1] += coeff_1;
  FdSFTranspose_vol_sym[2] += coeff_1;
  return 0;
};

static inline int ComputeFdSFTranspose_NeoHookean(double J_dVdJ, double J2_d2VdJ2, double lambda, double mu,
                                                               const double grad_du[3][3], double FdSFTranspose[3][3]) {
  double depsilon_sym[6], FdSFTranspose_vol_sym[6], FdSFTranspose_sym[6];
  // Compute depsilon = (grad_du + grad_du^T)/2
  depsilon_sym[0] = grad_du[0][0];
  depsilon_sym[1] = grad_du[1][1];
  depsilon_sym[2] = grad_du[2][2];
  depsilon_sym[3] = (grad_du[1][2] + grad_du[2][1]) / 2.;
  depsilon_sym[4] = (grad_du[0][2] + grad_du[2][0]) / 2.;
  depsilon_sym[5] = (grad_du[0][1] + grad_du[1][0]) / 2.;

  const double trace_depsilon = MatTraceSymmetric(depsilon_sym);

  // F*dS_vol*F^T
  ComputeFdSFTransposeVolumetric(J_dVdJ, J2_d2VdJ2, lambda, trace_depsilon, depsilon_sym, FdSFTranspose_vol_sym);

  // F*dS*F^T = F*dS_vol*F^T + 2 mu depsilon
  MatMatAddSymmetric(1.0, FdSFTranspose_vol_sym, 2 * mu, depsilon_sym, FdSFTranspose_sym);
  SymmetricMatUnpack(FdSFTranspose_sym, FdSFTranspose);
  return 0;
};

#ifdef __cplusplus
}
#endif

#endif // UTILS_H
