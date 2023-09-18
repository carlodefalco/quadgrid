#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

struct
compare {

  int count;
  compare () : count{0} {};

  bool
  operator()(int a, int b) {
    ++count;
    return a == b;
  }
  
};

struct
swap {

  int count;
  swap () : count{0} {};

  void
  operator()(int &a, int &b) {
    ++count;
    std::swap (a, b);
  }
  
};
  
int
main (int argc, char *argv[]) {

  compare c;
  swap s;
  
  constexpr int numelements = 1000;
  std::vector<int> start(numelements);
  for (int ii = 0; ii < numelements; ++ii) {
    start[ii] = ii;
  }
  std::vector<int> target = start;
  std::vector<int> ordering = start;

  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle (ordering.begin (), ordering.end (), g);

  std::cout << "ordering " << std::endl;
  for (auto const & kk : ordering) {
    std::cout << kk << std::endl;
  }
  
  for (int ii = 0; ii < target.size () - 1; ++ii) {

    if (! c (ii, ordering[ii])) {
      s (target[ii], target[ordering[ii]]);
      for (int jj = ii; jj < target.size (); ++jj) {
	if (c (ordering[jj], ii)) {
	  ordering[jj] = ordering[ii];
	  ordering[ii] = ii;
	  break;
	}
      }
    }
   
  }

  std::cout << "target " << std::endl;
  for (auto const & kk : target) {
    std::cout << kk << std::endl;
  }
  
  std::cout << "number of comparisons = " << c.count << std::endl;
  std::cout << "number of swaps = " << s.count << std::endl;
  return 0;
}

