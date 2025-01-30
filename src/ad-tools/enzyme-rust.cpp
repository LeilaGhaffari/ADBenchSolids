#include "../../include/ad-tools/enzyme-rust.h"

void init_data_enzyme_rust(double **stored_values, int Q) {
    *stored_values = (double *)malloc(Q * NUM_COMPONENTS_STORED_ENZYME_RUST * sizeof(double));
}

void free_data_enzyme_rust(double **stored_values) {
    if (*stored_values != NULL) {
        free(*stored_values);
        *stored_values = NULL;
    }
}

void tau_sym_enzyme_rust(const double e_sym[6], const double lambda, const double mu, double tau_sym[6]) {
    compute_stress(lambda, mu, e_sym, tau_sym);
}

// Residual Evaluation
void f_enzyme_rust(int Q, const double mu, const double lambda, double *dXdx_init, double *dudX, double **stored_values, double *f1) {
    double Grad_u[3][3], F_inv[3][3], tau_sym[6], tau[3][3], dXdx[3][3], e_sym[6];
    for (int i=0; i<Q; i++) {
        // Pack input data
        double dudX_loc[3][3], dXdx_init_loc[3][3];
        QDataPackMat(i, dXdx_init, dXdx_init_loc);
        QDataPackMat(i, dudX, dudX_loc);

        MatMatMult(1.0, dudX_loc, dXdx_init_loc, Grad_u);
        double F[3][3];
        DeformationGradient(Grad_u, F);
        const double Jm1  = MatDetAM1(Grad_u);
        const double detF = Jm1 + 1.0;
        MatInverse(F, detF, F_inv);
        MatMatMult(1.0, dXdx_init_loc, F_inv, dXdx);
        GreenEulerStrain(Grad_u, e_sym);
        tau_sym_enzyme_rust(e_sym, lambda, mu, tau_sym);
        SymmetricMatUnpack(tau_sym, tau);
        QDataUnpackMat(i, tau, f1);
        // Store
        StoredValuesPack(0, 9, NUM_COMPONENTS_STORED_ENZYME_RUST, i, (double *)dXdx, stored_values);
        StoredValuesPack(9, 6, NUM_COMPONENTS_STORED_ENZYME_RUST, i, (double *)e_sym, stored_values);
    }
}

void df_enzyme_rust(int Q, const double mu, const double lambda, double *ddudX, double **stored_values, double *df) {
}
