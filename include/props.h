//! \brief Class to represent a set of properties attached to a list of particles.

//! Data is stored in a contiguous memory region, each property
//! is a column, each particle a row. Extra space is (optionally)
//! allocated to allow for adding new particles/properties.
//! When a row is erased memory is not freed but kept for any subsequent
//! addition.

template <typename T>
struct
props_t {

  std::unique_ptr<T[]> buffer;
  
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
  };
  
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

  auto
  at (std::size_t iname) {   
    tcb::span retval (buffer.get () + iname * max_ptcls,
		      buffer.get () + iname * max_ptcls + num_ptcls);
    return retval;
  }
  
  void 
  erase_ptcl (std::size_t ir) {
    if (ir < this->num_ptcls) {
      for (auto const &prop : proplist) {
	auto col = this->at (prop.first);
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

