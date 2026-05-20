#include <tbasisfun.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace bspline;

int main(int argc, char *argv[]) {

  const char *out_path = (argc > 1) ? argv[1] : "tbasisfun_output.dat";
  std::ofstream out(out_path, std::ofstream::out);
  if (!out) {
    std::cerr << "Failed to open output file: " << out_path << std::endl;
    return 1;
  }

  int p{4};
  constexpr auto N = 30;
  int reg = 3;
  std::vector<double> b{0, 1, 2, 3, 4};
  std::vector<double> k = open_knot_vector (b.begin (), b.end (), p, reg);
  for (auto const & ii : k) std::cout << ii << std::endl;
  
  auto kb = k.begin ();
  auto ke = std::next (kb, p + 2);
  
  while (true) {
    
    std::vector<double> x(N+1, 0.0);
    std::vector<double> y(N+1, 0.0);

    x[0] = *kb;
    y[0] = onebasisfun(x[0], p, kb, ke);
    double dx = (*(std::prev (ke)) - *kb) / double (N);
    for (int ii = 1; ii < N+1; ++ii) {
      x[ii] = x[ii-1] + dx; 
      if (x[ii] >= (*(std::prev (ke)) - *(std::prev (ke)) * 4.0 * std::numeric_limits<double>::epsilon ()))
	x[ii] = *(std::prev (ke));
      y[ii] = onebasisfun (x[ii], p, kb, ke);
    }

    for (int ii = 0; ii < N+1; ++ii) {
      out << std::setprecision(19) << x[ii] << " " << y[ii] << "\n";
    }
    out << "\n";
    
    if (ke != k.end ()) {
      ++kb; ++ke;
    }
    else {
      break;
    }
  }
  
  //std::cout  << std::setprecision (19) << *std::prev(ke) << "   " << onebasisfun<Position::Boundary> (*std::prev(ke) , p, kb, ke) << std::endl;

  out.close();

  return 0;  
}

