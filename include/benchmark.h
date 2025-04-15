#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "get-data.h"

void DisplayTimeAndError(const std::vector<std::string> &ad_tools, int Q,
                         std::vector<double> &dXdx_init,
                         std::vector<double> &dudX,
                         std::vector<double> &ddudX) {
  // Header
  int tool_width = 20, time_width = 21, error_width = 20, line_width = 91;
  std::cout << std::string(line_width, '-') << std::endl;
  std::cout << std::left << std::setw(tool_width) << "AD Tool"
            << std::setw(time_width) << "Residual Time (s)"
            << std::setw(error_width) << "Residual Error"
            << std::setw(time_width) << "Jacobian Time (s)"
            << std::setw(error_width) << "Jacobian Error" << std::endl;
  std::cout << std::string(line_width, '-') << std::endl;

  // analytic solution
  const double mu = 1., lambda = 1.0;
  Bench bench_ref;
  const std::string &analytic = "analytic-c";
  bench_setup(&bench_ref, analytic.c_str());
  double *stored_values_ref = NULL;
  int n;
  bench_ref.init_data(&stored_values_ref, Q, &n);
  double *f_ref = (double *)malloc(Q * 9 * sizeof(double));
  double *df_ref = (double *)malloc(Q * 9 * sizeof(double));
  bench_ref.f(Q, mu, lambda, dXdx_init.data(), dudX.data(), &stored_values_ref,
              f_ref);
  bench_ref.df(Q, mu, lambda, ddudX.data(), &stored_values_ref, df_ref);

  // Loop over models/AD tools
  for (const auto &ad_tool : ad_tools) {
    Bench bench;
    if (bench_setup(&bench, ad_tool.c_str()) != 0) {
      std::cerr << "Failed to set up bench for AD tool: " << ad_tool
                << std::endl;
      return;
    }
    double *stored_values = NULL;
    double *f = (double *)malloc(Q * 9 * sizeof(double));
    double *df = (double *)malloc(Q * 9 * sizeof(double));
    double f_total_error = 0.0, df_total_error = 0.0;
    int num_comp;
    bench.init_data(&stored_values, Q, &num_comp);

    // Measure time for f
    // Run before timing to make sure the data is hot
    bench.f(Q, mu, lambda, dXdx_init.data(), dudX.data(), &stored_values, f);
    auto start_f = std::chrono::high_resolution_clock::now();
    bench.f(Q, mu, lambda, dXdx_init.data(), dudX.data(), &stored_values, f);
    auto end_f = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_f = end_f - start_f;

    // Measure error for f
    f_total_error += ComputeL2Error(f, f_ref, Q * 9);

    // Measure time for df
    std::chrono::duration<double> elapsed_df{};

    // Run before timing to make sure the data is hot
    bench.df(Q, mu, lambda, ddudX.data(), &stored_values, df);
    auto start_df = std::chrono::high_resolution_clock::now();
    bench.df(Q, mu, lambda, ddudX.data(), &stored_values, df);
    auto end_df = std::chrono::high_resolution_clock::now();
    elapsed_df = end_df - start_df;

    // Measure error for df
    df_total_error += ComputeL2Error(df, df_ref, Q * 9);

    // Print results
    if (ad_tool == "stream-triad" || ad_tool == "stream-residual") {
      std::cout << std::left << std::setw(tool_width) << ad_tool
                << std::setw(time_width) << elapsed_f.count() << std::endl;
    } else {
      std::cout << std::left << std::setw(tool_width) << ad_tool
                << std::setw(time_width) << elapsed_f.count()
                << std::setw(error_width) << f_total_error
                << std::setw(time_width) << elapsed_df.count()
                << std::setw(error_width) << df_total_error << std::endl;
    }

    // Cleanup
    free(f);
    free(df);
    bench.free_data(&stored_values);
  }
  std::cout << std::string(line_width, '-') << std::endl;

  // Cleanup
  free(f_ref);
  free(df_ref);
  bench_ref.free_data(&stored_values_ref);
}

#endif // BENCHMARK_H
