#ifndef GET_DATA_H
#define GET_DATA_H

std::vector<std::string> parse_model_arg(const std::string &input) {
  std::vector<std::string> models;
  std::stringstream ss(input);
  std::string model;
  while (std::getline(ss, model, ',')) {
    if (!model.empty()) {
      models.push_back(model);
    }
  }
  return models;
}

void GetData(const std::string &filename, int Q, std::vector<double> &dXdx_init,
             std::vector<double> &dudX, std::vector<double> &ddudX) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Unable to open file\n";
    std::exit(EXIT_FAILURE);
  }

  std::string line;
  std::getline(file, line); // Skip the header row

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
    dudX[idx] = value2;
    ddudX[idx] = value3;

    ++idx;
  }
  file.close();
}

int QuadraturePointsNumber(const std::string &filename) {
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

// For Debugging
void PrintData(const double *f, int n) {
  std::cout << "\n";
  for (int i = 0; i < n; ++i) {
    std::cout << f[i] << std::endl;
  }
  std::cout << "\n";
}

// Transpose data from [N, 3, 3] to [3, 3, N]
void TransposeData(std::vector<double> &in, std::vector<double> &out, int Q) {
  for (int i = 0; i < Q; i++)
    for (int j = 0; j < 9; j++)
      out[i + Q * j] = in[i * 9 + j];
}

#endif // GET_DATA_H
