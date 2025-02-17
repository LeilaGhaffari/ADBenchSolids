#include "../../include/ad-tools/analytical.h"
#include <omp.h>

void init_data_analytic(double **stored_values, int Q) {
    *stored_values = (double *)malloc(Q * NUM_COMPONENTS_STORED_ANALYTICAL * sizeof(double));
}

void free_data_analytic(double **stored_values) {
    if (*stored_values != NULL) {
        free(*stored_values);
        *stored_values = NULL;
    }
}

void f_analytic(int Q, const double mu, const double lambda, double *dXdx_init, double *dudX,
                double **stored_values, double *f1) {
    #pragma omp parallel for
    for (int i=0; i<Q; i++) {
        double Grad_u[3][3], F_inv[3][3], tau[3][3], tau_sym[6], dXdx[3][3], e_sym[6],
               J_dVdJ, F[3][3], dudX_loc[3][3], dXdx_init_loc[3][3];
        // Pack input data
        QDataPackMat(i, Q, dXdx_init, dXdx_init_loc);
        QDataPackMat(i, Q, dudX, dudX_loc);

        // Grad_u = du/dx_initial = du/dX * dX/dx_initial
        MatMatMult(1.0, dudX_loc, dXdx_init_loc, Grad_u);

        // Compute the Deformation Gradient : F = I + Grad_u
        DeformationGradient(Grad_u, F);
        GreenEulerStrain(Grad_u, e_sym);
        const double Jm1 = MatDetAM1(Grad_u);
        VolumetricFunctionAndDerivatives(Jm1, NULL, &J_dVdJ, NULL);
        KirchhoffTau_NeoHookean(J_dVdJ, lambda, 2*mu, e_sym, tau_sym);
        SymmetricMatUnpack(tau_sym, tau);
        QDataUnpackMat(i, Q, tau, f1);

        // Compute F^{-1}
        const double detF = Jm1 + 1.;
        MatInverse(F, detF, F_inv);

        // dXdx = dX/dx = dX/dx_initial * F^{-1}
        MatMatMult(1.0, dXdx_init_loc, F_inv, dXdx);

        // Store
        StoredValuesPack(Q, i, 0, 9, (double *)dXdx, stored_values);
        StoredValuesPack(Q, i, 9, 6, (double *)tau_sym, stored_values);
        StoredValuesPack(Q, i, 15, 1, &Jm1, stored_values);
    }
}

void df_analytic(int Q, const double mu, const double lambda,
                 double *ddudX, double **stored_values, double *df) {
    #pragma omp parallel for
    for (int i=0; i<Q; i++) {
        double grad_du[3][3], tau_sym[6], tau[3][3], grad_du_tau[3][3], df_mat[3][3],
               FdSFTranspose[3][3], dXdx[3][3], Jm1, J_dVdJ, J2_d2VdJ2, ddudX_loc[3][3];
        // Unpack stored values
        StoredValuesUnpack(Q, i, 0, 9, stored_values, (double *)dXdx);
        StoredValuesUnpack(Q, i, 9, 6, stored_values, (double *)tau_sym);
        StoredValuesUnpack(Q, i, 15, 1, stored_values, &Jm1);

        // Pack input data
        QDataPackMat(i, Q, ddudX, ddudX_loc);

        // Compute grad_du = ddu/dX * dX/dx
        MatMatMult(1.0, ddudX_loc, dXdx, grad_du);
        SymmetricMatUnpack(tau_sym, tau);

        // Compute grad_du_tau = grad_du*tau
        MatMatMult(1.0, grad_du, tau, grad_du_tau);

        // Compute F*dS*F^T
        VolumetricFunctionAndDerivatives(Jm1, NULL, &J_dVdJ, &J2_d2VdJ2);
        ComputeFdSFTranspose_NeoHookean(J_dVdJ, J2_d2VdJ2, lambda, mu, grad_du, FdSFTranspose);

        // df1 = dtau - tau * grad_du^T
        //     = grad_du*tau + F*dS*F^T
        MatMatAdd(1.0, grad_du_tau, 1., FdSFTranspose, df_mat);
        QDataUnpackMat(i, Q, df_mat, df);
    }
}
