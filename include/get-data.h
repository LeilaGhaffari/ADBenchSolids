#ifndef GET_DATA_H
#define GET_DATA_H

// For Debugging
void PrintData(const double *f, int n) {
    std::cout << "\n";
    for (int i = 0; i < n; ++i) {
        std::cout << f[i] << std::endl;
    }
    std::cout << "\n";
}

// Transpose data from [N, 3, 3] to [3, 3, N]
void TransposeData(std::vector<double>& in, std::vector<double>& out, int Q) {
    for (int i=0; i<Q; i++) for (int j=0; j<9; j++) out[i + Q*j] = in[i*9 + j];
}

void GetData(const std::string& filename, int Q, std::vector<double>& dXdx_init,
             std::vector<double>& dudX, std::vector<double>& ddudX) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file\n";
        std::exit(EXIT_FAILURE);
    }

    std::string line;
    std::getline(file, line);  // Skip the header row

    int idx = 0;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        double value1, value2, value3;
        std::getline(ss, value, ',');
        value1 = std::stod(value);
        std::getline(ss, value, ',');
        value2 = std::stod(value);
        std::getline(ss, value, ',');
        value3 = std::stod(value);

        dXdx_init[idx] = value1;
        dudX[idx]= value2;
        ddudX[idx] = value3;

        ++idx;
    }
    file.close();
}

int QuadraturePointsNumber(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::string line;
    int row_count = 0;

    std::getline(file, line); // Skip the header
    while (std::getline(file, line)) {
        ++row_count;
    }

    file.close();
    return row_count / 9;
}

double ComputeL2Error(const double *f, const double *f_ref, int n) {
    double error = 0.0;
    for (int i = 0; i < n; ++i) {
        double diff = f[i] - f_ref[i];
        error += diff * diff;
    }
    return std::sqrt(error);
}

void DisplayTimeAndError(const std::vector<std::string> &ad_tools,
                             int Q,
                             std::vector<double>& dXdx_init,
                             std::vector<double>& dudX,
                             std::vector<double>& ddudX) {
    // Header
    int tool_width = 15, time_width = 21, error_width = 20, line_width = 91;
    std::cout << std::string(line_width, '-') << std::endl;
    std::cout << std::left
              << std::setw(tool_width) << "AD Tool"
              << std::setw(time_width) << "Residual Time (s)"
              << std::setw(error_width) << "Residual Error"
              << std::setw(time_width) << "Jacobian Time (s)"
              << std::setw(error_width) << "Jacobian Error"
              << std::endl;
    std::cout << std::string(line_width, '-') << std::endl;

    // Analytical solution
    const double mu = 1., lambda = 1.0;
    Bench bench_ref;
    const std::string &analytical = "analytical";
    bench_setup(&bench_ref, analytical.c_str());
    double *stored_values_ref = NULL;
    bench_ref.init_data(&stored_values_ref, Q);
    double *f_ref = (double *)calloc(Q * 9, sizeof(double));
    double *df_ref = (double *)calloc(Q * 9, sizeof(double));
    bench_ref.f(Q, mu, lambda, dXdx_init.data(), dudX.data(), &stored_values_ref, f_ref);
    bench_ref.df(Q, mu, lambda, ddudX.data(), &stored_values_ref, df_ref);

    // Loop over models/AD tools
    for (const auto &ad_tool : ad_tools) {
        Bench bench;
        if (bench_setup(&bench, ad_tool.c_str()) != 0) {
            std::cerr << "Failed to set up bench for AD tool: " << ad_tool << std::endl;
            return;
        }
        double *stored_values = NULL;
        double *f = (double *)calloc(Q * 9, sizeof(double));
        double *df = (double *)calloc(Q * 9, sizeof(double));
        double f_total_error = 0.0, df_total_error = 0.0;
        bench.init_data(&stored_values, Q);

        // Measure time for f
        auto start_f = std::chrono::high_resolution_clock::now();
        bench.f(Q, mu, lambda, dXdx_init.data(), dudX.data(), &stored_values, f);
        auto end_f = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_f = end_f - start_f;

        // Measure error for f
        f_total_error += ComputeL2Error(f, f_ref, Q*9);

        // Measure time for df
        std::chrono::duration<double> elapsed_df{};
        if (ad_tool == "enzyme-rust") {
            df_total_error = std::numeric_limits<double>::infinity();
        } else {
            auto start_df = std::chrono::high_resolution_clock::now();
            bench.df(Q, mu, lambda, ddudX.data(), &stored_values, df);
            auto end_df = std::chrono::high_resolution_clock::now();
            elapsed_df = end_df - start_df;

            // Measure error for df
            df_total_error += ComputeL2Error(df, df_ref, Q*9);
        }

        // Print results
        std::cout << std::left << std::setw(tool_width) << ad_tool
                  << std::setw(time_width) << elapsed_f.count()
                  << std::setw(error_width) << f_total_error
                  << std::setw(time_width) << elapsed_df.count()
                  << std::setw(error_width) << df_total_error
                  << std::endl;

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

#endif // GET_DATA_H
