#include <tbasisfun.h>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace bspline;

int main(int argc, char *argv[]) {

  int p{1};
  constexpr auto N = 30;
  
  std::vector<double> b{0.0, 1.0, 2.5, 4.0, 7.0};
  std::vector<double> k = open_knot_vector (b.begin (), b.end (), p, p - 1);
  
  auto kb = k.begin ();
  auto ke = std::next (kb, p + 2);
  
  while (true) {
    
    std::vector<double> x(N+1, 0.0);
    std::vector<double> y(N+1, 0.0);

    x[0] = *kb;
    y[0] = onebasisfun (x[0], p, kb, ke);
    double dx = (*(std::prev (ke)) - *kb) / double (N);
    for (int ii = 1; ii < N+1; ++ii) {
      x[ii] = x[ii-1] + dx; 
      if (x[ii] >= (*(std::prev (ke)) - *(std::prev (ke)) * 4.0 * std::numeric_limits<double>::epsilon ()))
	x[ii] = *(std::prev (ke));
      y[ii] = onebasisfunder (x[ii], p, kb, ke);
    }

    for (int ii = 0; ii < N+1; ++ii) {
      std::cout  << std::setprecision (15) << x[ii] << "   " << y[ii] << std::endl;
    }
    
    if (ke != k.end ()) {
      ++kb;
      ++ke;
    }
    else {
      break;
    }
  }
  
  return 0;  
}

