#include <tbasisfun.h>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace bspline;

int main(int argc, char *argv[]) {

  int p{2};
  constexpr auto N = 100;
  
  std::vector<double> k{1.0, 1.99, 2.01, 3.0};
  auto kb = k.begin ();
  auto ke = k.end ();

  /*
  std::vector<double> x(N+1, 0.0);
  std::vector<double> y(N+1, 0.0);

  x[0] = *kb;
  y[0] = onebasisfun (x[0], p, kb, ke);
  double dx = (*(std::prev (ke)) - *kb) / double (N);
  for (int ii = 1; ii < N+1; ++ii) {
    x[ii] = x[ii-1] + dx; 
    if (x[ii] >= (*(std::prev (ke)) - *(std::prev (ke)) * 4.0 * std::numeric_limits<double>::epsilon ()))
      x[ii] = *(std::prev (ke));
    y[ii] = onebasisfun (x[ii], p, kb, ke);
  }

  for (int ii = 0; ii < N+1; ++ii) {
    std::cout  << std::setprecision (19) << x[ii] << "   " << y[ii] << std::endl;
  }
  */

  double x{2.0};
  double y  = onebasisfun (x, p, kb, ke);
  std::cout  << std::setprecision (19) << x << "   " << y << std::endl;
  
  return 0;  
}

