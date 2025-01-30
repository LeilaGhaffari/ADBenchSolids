#ifndef BENCH_H
#define BENCH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <cstring>
#include "ad-tools/adolc.h"
#include "ad-tools/enzyme.h"
#include "ad-tools/enzyme-rust.h"
#include "ad-tools/analytical.h"
#include "ad-tools/tapenade.h"

#define SETUP_BENCH(tool_name)                    \
    bench->init_data = init_data_##tool_name;     \
    bench->free_data = free_data_##tool_name;     \
    bench->f = f_##tool_name;                     \
    bench->df = df_##tool_name;

typedef struct Bench {
    void (*init_data)(double **stored_values, int Q);
    void (*free_data)(double **stored_values);
    void (*f)(int Q, const double mu, const double lambda, double *dXdx_init, double *dudX, double **stored_values, double *f1);
    void (*df)(int Q, const double mu, const double lambda, double *ddudX, double **stored_values, double *df);
} Bench;

int bench_setup(Bench *bench, const char *tool) {
    if (strcmp(tool, "adolc") == 0) {
        SETUP_BENCH(adolc);
    } else if (strcmp(tool, "enzyme-c") == 0) {
        SETUP_BENCH(enzyme);
    } else if (strcmp(tool, "enzyme-rust") == 0) {
        SETUP_BENCH(enzyme_rust);
    } else if (strcmp(tool, "tapenade") == 0) {
        SETUP_BENCH(tapenade);
    } else if (strcmp(tool, "analytical") == 0) {
        SETUP_BENCH(analytic);
    } else {
        printf("Unknown model: %s\n", tool);
        printf("Valid options are: analytical, adolc, enzyme-c, enzyme-rust, and tapenade\n");
        return 1;
    }
    return 0;
}

#endif // BENCH_H
