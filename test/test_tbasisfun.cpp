#include <tbasisfun.h>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace bspline;

int main(int argc, char *argv[]) {

  int p{3};
  constexpr auto N = 1000;
  std::vector<double> x(N+1, 0.0);
  std::vector<double> y(N+1, 0.0);
  std::vector<double> U{0.0, 0.0, 4.5, 9.0, 9.0};

  x[0] = -1.;
  double dx = 12. / double (N);
  for (int ii = 1; ii < N+1; ++ii) {
    x[ii] = dx * ii;
    y[ii] = onebasisfun (x[ii], p, U.cbegin (), U.cend ());
  }
  
  for (int ii = 0; ii < N+1; ++ii) {
    std::cout << std::setprecision (16) << x[ii] << "   " << y[ii] << std::endl;
  }
 
  return 0;  
}

