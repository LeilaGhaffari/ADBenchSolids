#include "../../include/ad-tools/stream-residual.h"

void init_data_stream_residual(double **stored_values, int Q, int *num_comp) {
  *stored_values = (double *)malloc(Q * NUM_COMPONENTS_STORED_STREAM_RESIDUAL *
                                    sizeof(double));
  *num_comp = NUM_COMPONENTS_STORED_STREAM_RESIDUAL;
}

void free_data_stream_residual(double **stored_values) {
  if (*stored_values != NULL) {
    free(*stored_values);
    *stored_values = NULL;
  }
}

void f_stream_residual(int Q, const double mu, const double lambda,
                       double *dXdx_init, double *dudX, double **stored_values,
                       double *f1) {
  BenchPragmaSIMD for (int i = 0; i < Q; i++) {
    double A[3][3], B[3][3], C[3][3];
    // Pack input data
    QDataPackMat(i, Q, dXdx_init, A);
    QDataPackMat(i, Q, dudX, B);

    // C = mu * A + B
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        C[j][k] = mu * A[j][k] + B[j][k];
      }
    }
    QDataUnpackMat(i, Q, C, f1);

    // Store (to make it comparable to analytic-c)
    StoredValuesPack(Q, i, 0, 9, (double *)A, stored_values);
    StoredValuesPack(Q, i, 9, 7, (double *)B, stored_values);
  }
}

void df_stream_residual(int Q, const double mu, const double lambda,
                        double *ddudX, double **stored_values, double *df) {}
