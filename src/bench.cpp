#include "../include/bench.h"
#include "../include/get-data.h"

int main(int argc, char *argv[]) {
  // AD tools
  std::vector<std::string> ad_tools = {"analytical", "enzyme-c", "enzyme-rust", "tapenade", "adolc"};

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

  // Compute Time and Error in evaluating residual and Jacobian
  DisplayTimeAndError(ad_tools, Q, dXdx_init, dudX, ddudX);

  return 0;
}
