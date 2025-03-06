#ifndef BENCH_H
#define BENCH_H

#include "ad-tools/adolc.h"
#include "ad-tools/analytic-c.h"
#include "ad-tools/analytic-rust.h"
#include "ad-tools/enzyme-c.h"
#include "ad-tools/enzyme-rust.h"
#include "ad-tools/tapenade.h"
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#define SETUP_BENCH(tool_name)                                                 \
  bench->init_data = init_data_##tool_name;                                    \
  bench->free_data = free_data_##tool_name;                                    \
  bench->f = f_##tool_name;                                                    \
  bench->df = df_##tool_name;

typedef struct Bench {
  void (*init_data)(double **stored_values, int Q, int *num_comp);
  void (*free_data)(double **stored_values);
  void (*f)(int Q, const double mu, const double lambda, double *dXdx_init,
            double *dudX, double **stored_values, double *f1);
  void (*df)(int Q, const double mu, const double lambda, double *ddudX,
             double **stored_values, double *df);
} Bench;

int bench_setup(Bench *bench, const char *tool) {
  if (strcmp(tool, "analytic-c") == 0) {
    SETUP_BENCH(analytic_c);
  } else if (strcmp(tool, "analytic-rust") == 0) {
    SETUP_BENCH(analytic_rust);
  } else if (strcmp(tool, "enzyme-c") == 0) {
    SETUP_BENCH(enzyme_c);
  } else if (strcmp(tool, "enzyme-rust") == 0) {
    SETUP_BENCH(enzyme_rust);
  } else if (strcmp(tool, "tapenade") == 0) {
    SETUP_BENCH(tapenade);
  } else if (strcmp(tool, "adolc") == 0) {
    SETUP_BENCH(adolc);
  } else {
    printf("Unknown model: %s\n", tool);
    printf("Valid options are: analytic-c, analytic-rust, enzyme-c, "
           "enzyme-rust, tapenade, and adolc\n");
    return 1;
  }
  return 0;
}

#endif // BENCH_H
