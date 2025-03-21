#include "../../include/ad-tools/stream.h"

void init_data_stream(double **stored_values, int Q, int *num_comp) {
  *stored_values =
      (double *)malloc(Q * NUM_COMPONENTS_STORED_STREAM * sizeof(double));
  *num_comp = NUM_COMPONENTS_STORED_STREAM;
}

void free_data_stream(double **stored_values) {
  if (*stored_values != NULL) {
    free(*stored_values);
    *stored_values = NULL;
  }
}

void f_stream(int Q, const double mu, const double lambda, double *dXdx_init,
              double *dudX, double **stored_values, double *f1) {
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
  }
}

void df_stream(int Q, const double mu, const double lambda, double *ddudX,
               double **stored_values, double *df) {}
