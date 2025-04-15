#include "../include/bench.h"
#include "../include/benchmark.h"
#include "../include/get-data.h"

int main(int argc, char *argv[]) {
  // Process CLI
  std::vector<std::string> ad_tools = {
      "stream-triad", "stream-residual", "analytic-c", "analytic-rust",
      "enzyme-c",     "enzyme-rust",     "tapenade",   "adolc"};
  std::string filename = "";
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg.rfind("-models=", 0) == 0) {
      std::string list = arg.substr(8);
      ad_tools = parse_model_arg(list);
    } else if (arg.rfind("-data=", 0) == 0) {
      filename = arg.substr(6);
    }
  }

  // Get Q
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

  // Free memory used by old vectors
  std::vector<double>().swap(dXdx_init);
  std::vector<double>().swap(dudX);
  std::vector<double>().swap(ddudX);

  // Compute Time and Error in evaluating residual and Jacobian
  DisplayTimeAndError(ad_tools, Q, dXdx_initT, dudXT, ddudXT);

  return 0;
}
