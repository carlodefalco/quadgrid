#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <tcb/span.hpp>
#include <vector>

template <typename T>
struct
props_t {

  std::unique_ptr<T[]> buffer;
  using difference_type = decltype (buffer.get () - buffer.get ());
  
  const std::size_t max_props; /*!< maximun number of T-valued properties */
  std::size_t num_props;       /*!< actual current number of T-valued properties */

  const std::size_t max_ptcls; /*!< maximun number of particles */
  std::size_t num_ptcls;       /*!< actual current numberof particles */

  std::map<std::string, std::size_t> proplist{};
  
  props_t () = delete;
  props_t (const props_t&) = delete;

  void
  init_names (const std::vector<std::string>& propnames) {
    std::size_t ii = 0;
    for (auto & j : propnames)
      proplist[j] = ii++;
  }
  
  props_t (std::size_t max_props_, std::size_t num_props_,
	    std::size_t max_ptcls_, std::size_t num_ptcls_)
    :  max_props(max_props_), num_props(num_props_),
       max_ptcls(max_ptcls_), num_ptcls(num_ptcls_) {
    buffer = std::make_unique<T[]> (max_ptcls * max_props);    
  }
  

  props_t (std::size_t max_props_, const std::vector<std::string>& propnames,
	    std::size_t max_ptcls_, std::size_t num_ptcls_)
    : props_t(max_props_, propnames.size (), max_ptcls_, num_ptcls_) {
    init_names (propnames);    
  }

  props_t (const std::vector<std::string>& propnames, std::size_t max_ptcls_)
    : props_t (propnames.size (), propnames, max_ptcls_, max_ptcls_) { }

  
  auto
  at (const std::string &name) {   
    tcb::span retval (buffer.get () + proplist.at (name) * max_ptcls,
		      buffer.get () + proplist.at (name) * max_ptcls + num_ptcls);
    return retval;
  }

  void 
  erase_ptcl (std::size_t ir) {
    tcb::span<T> col;
    if (ir < this->num_ptcls) {
      for (auto const &prop : proplist) {
	col = this->at (prop.first);
	std::copy (std::next (col.begin (), ir + 1), col.end (), std::next (col.begin (), ir));
      }
      --(this->num_ptcls);
    } else {
      throw std::out_of_range ("particle index too large");
    }    
  }

  
  void
  reorder (std::vector<std::size_t> ordering) {
    tcb::span<T> col;
   
    for (std::size_t ii = 0; ii < num_ptcls - 1; ++ii) {

      if (ii != ordering[ii]) {
	for (auto const &prop : proplist) {
	  col = this->at (prop.first);
	  std::swap (col[ii], col[ordering[ii]]);
	}
	for (int jj = ii; jj < col.size (); ++jj) {
	  if (ordering[jj] == ii) {
	    ordering[jj] = ordering[ii];
	    ordering[ii] = ii;
	    break;
	  }
	}
      }
      
    }
  }
  
};

using dprops_t = props_t<double>;
using iprops_t = props_t<int>;

int
main () {
  
  std::vector<std::string> list{"a", "b", "c"};
  dprops_t d (list, 10);

  auto a = d.at ("a");
  std::cout << "size (a) = " << a.size () << std::endl;
  
  double ii{0.0};
  for (auto & ia : a) {
    ia = ii;
    ii += 1.1;
  }

  for (auto const & ia : a) {
    std::cout << ia << std::endl;
  }

  d.erase_ptcl (d.num_ptcls - 2);
  a = d.at ("a");
  std::cout << "size (a) = " << a.size () << std::endl;
  std::cout << "a = " << std::endl;
  for (auto const & ia : a) {
    std::cout << ia << std::endl;
  }

  std::vector<std::size_t> start (d.num_ptcls);
  for (int ii = 0; ii < start.size (); ++ii) {
    start[ii] = ii;
  }
  std::vector<std::size_t> ordering = start;

  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle (ordering.begin (), ordering.end (), g);

  std::cout << "ordering = " << std::endl;
  for (auto const & ia : ordering) {
    std::cout << ia << std::endl;
  }

  d.reorder (ordering);
  a = d.at ("a");
  std::cout << "a = " << std::endl;
  for (auto const & ia : a) {
    std::cout << ia << std::endl;
  }
  
  return 0;
}

