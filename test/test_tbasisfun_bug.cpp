#include <tbasisfun.h>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace bspline;

int main(int argc, char *argv[]) {

  int p{2};

  // This is a knot vector with a knot of multiplicity p
  // so that ath the corresponding break the funcion is C0
  // and interpolating.
  std::vector<double> k{1.0, 2.0, 2.0, 3.0};
  
  auto kb = k.begin ();
  auto ke = k.end ();

  // Evaluate the basis function EXACTLY at the break point
  double x{2.0};
  double y  = onebasisfun (x, p, kb, ke);

  // The value of the basis function should be 1.0, but is 2.0 if the
  // code does not understand that there is no knot span between two
  // identical knots.
  std::cout  << std::setprecision (19) << "x = (should be 2.0) " << x << " y = (should be 1.0)  " << y << std::endl;
  
  return 0;  
}

