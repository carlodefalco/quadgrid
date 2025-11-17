#include <mmorton.h>
#include <array>
#include <vector>
#include <iostream>

int main ()
{

  std::cout << "using " << sizeof(morton_code_t)
	    << " bytes for morton numbers" << std::endl;

  std::cout << "using " << sizeof(coord_t)
	    << " bytes for coordinates" << std::endl;

  uint ndivs = 32;
  std::vector<std::array<int, 2>> v(ndivs*ndivs);
  for (int ii = 0; ii < ndivs; ++ii)
    for (int jj = 0; jj < ndivs; ++jj) {
      auto kk = coord_2_morton (ii, jj);
      v[kk] = {ii, jj};
    }

  // Human readable
  /*  for (auto dd = 0; dd < v.size (); ++dd) {
      std::cout << "(x, y) = (" << v[dd][0] << ", " << v[dd][1]
      << "), mn = " << dd;
      coord_t x = 0, y = 0;
      morton_2_coord (dd, x, y);
      std::cout << " (x, y) = (" << x << ", " << y
      << ")" << std::endl;
      }
  */

  // File for plotting
  for (auto dd = 0; dd < v.size (); ++dd) {
    std::cout << dd << ", " << v[dd][0]
	      << ", " << v[dd][1]
	      <<  std::endl;
  }
  
  /*std::cout << coord_2_morton (3, 0) << std::endl;*/
  return 0;
  
};
