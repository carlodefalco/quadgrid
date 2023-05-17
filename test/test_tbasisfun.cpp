#include <tbasisfun.h>
#include <iomanip>
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {

  int p{2};
  std::vector<double> x(101, 0.0);
  std::vector<double> y(101, 0.0);
  std::vector<double> U{0.0, 1.0, 5.0, 9.0};

  x[0] = U[0];
  double dx = (U[3] - U[0]) / 100.0;
  for (int ii = 1; ii < 101; ++ii) {
    x[ii] = dx * ii;
    y[ii] = onebasisfun (x[ii], 2, U.cbegin (), U.cend ());
  }

  for (int ii = 0; ii < 101; ++ii) {
    std::cout << std::setprecision (16) << x[ii] << "   " << y[ii] << std::endl;
  }
  return 0;  
}

