#include "tapenade.h"
#include <cstdio>
#include <iostream>

int main() {
  TapenadeContext *context = new TapenadeContext;
  context->mu = 1.;
  context->lambda = 1.;
  context->stored = new double[NUM_COMPONENTS_STORED_TAPENADE];

  // Residual Evaluation
  const double dXdx_initial[3][3] = {{0.0702417, 0.4799115, 0.3991242},
                                     {0.6756593, 0.0633284, 0.0959267},
                                     {0.2241923, 0.0281781, 0.0917613}};
  const double dudX[3][3] = {{0.1425560, 0.115120, 0.551640},
                             {0.0591922, 0.123535, 0.166572},
                             {0.1617210, 0.478828, 0.646217}};
  double f1[3][3];
  f_tapenade(context, dXdx_initial, dudX, f1);

  // Jacobian Evaluation
  double ddudX[3][3] = {{0.23617975, 0.60250516, 0.1717169},
                        {0.86615524, 0.3365063, 0.17371375},
                        {0.0441905, 0.16762188, 0.45047923}};
  double df1[3][3];
  df_tapenade(context, ddudX, df1);

  free(context->stored);
  free(context);

  return 0;
}
