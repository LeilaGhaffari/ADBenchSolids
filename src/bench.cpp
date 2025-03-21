#include "../include/bench.h"
#include "../include/get-data.h"

std::vector<std::string> parse_model_arg(const std::string &input) {
  std::vector<std::string> models;
  std::stringstream ss(input);
  std::string model;
  while (std::getline(ss, model, ',')) {
    if (!model.empty()) {
      models.push_back(model);
    }
  }
  if (models.empty()) {
    models.push_back("stream");
  }
  return models;
}

int main(int argc, char *argv[]) {
  // AD tools
  std::vector<std::string> ad_tools = {"stream"}; // Default model
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg.rfind("-models=", 0) == 0) {
      std::string list = arg.substr(8); // Remove "-models="
      ad_tools = parse_model_arg(list);
    }
  }

  // File setup
  std::string filename = "random-data.csv";
  int mat_elem = 9, Q = QuadraturePointsNumber(filename);
  std::cout << "\nQuadrature Points = " << Q << "\n\n";

  // Allocate input arrays as 1D vectors
  std::vector<double> dXdx_init(mat_elem * Q, 0.0);
  std::vector<double> dudX(mat_elem * Q, 0.0);
  std::vector<double> ddudX(mat_elem * Q, 0.0);

  // Read the data from the CSV file
  GetData(filename, Q, dXdx_init, dudX, ddudX);

  // Transpose data
  std::vector<double> dXdx_initT(mat_elem * Q, 0.0);
  std::vector<double> dudXT(mat_elem * Q, 0.0);
  std::vector<double> ddudXT(mat_elem * Q, 0.0);
  TransposeData(dXdx_init, dXdx_initT, Q);
  TransposeData(dudX, dudXT, Q);
  TransposeData(ddudX, ddudXT, Q);

  // Compute Time and Error in evaluating residual and Jacobian
  DisplayTimeAndError(ad_tools, Q, dXdx_initT, dudXT, ddudXT);

  return 0;
}
