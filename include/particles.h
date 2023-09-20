#ifndef PARTICLES_H
#define PARTICLES_H

#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <json.hpp>
#include <map>
#include <quadgrid_cpp.h>
#include <string>

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
s  at (const std::string &name) {   
    tcb::span retval (buffer.get () + proplist.at (name) * max_ptcls,
		      buffer.get () + proplist.at (name) * max_ptcls + num_ptcls);
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

using dprops_t = props_t<double>;
using iprops_t = props_t<int>;


//! datatype for assignment operators
using assignment_t = std::function <double& (double&, const double&)>;

namespace ASSIGNMENT_OPS {
  extern assignment_t EQ; // = [] (double& TO, const double& FROM) -> double& { return TO = FROM; };
  extern assignment_t PLUS_EQ; // = [] (double& TO, const double& FROM) -> double& { return TO += FROM; };
  extern assignment_t TIMES_EQ; // = [] (double& TO, const double& FROM) -> double& { return TO *= FROM; };
}

//! \brief Class to represent particles embedded in a grid.

//! Offers methods for transfer of quantities from particles
//! to grid or vice-versa (particles_t::p2g, `g2p, p2gd`, `g2pd`).
//! Initial positions of the particles are chosen at random
//! unless otherwise specified.
//! To each particle a set one can associate a set of `double`
//! and one of `int` which are stored in the `std::map` variables
//! `dprops` and `iprops`, respectively.
//! Can compute a (lumped) mass matrix to be used in the transfer
//! functions.
//! If particles are moved, the connectivity must be updated invoking
//! the method `init_particle_mesh ()`

struct
particles_t {

  //! datatype for indexing into vectors of properties
  using idx_t = std::size_t;
  
  idx_t num_particles;    //!< number of particles.
  std::vector<double> x;  //!< x coordinate of particle positions.
  std::vector<double> y;  //!< y coordinate of particle positions.

  //! integer type quantities associated with the particles.
  iprops_t iprops;
  
  //! `double` type quantities associated with the particles.
  dprops_t dprops;  

  std::vector<double> M; //!< Mass matrix to be used for transfers if required.
  std::map<idx_t, std::vector<idx_t>> grd_to_ptcl;   //!< grid/particles connectivity.
  const quadgrid_t<std::vector<double>>& grid;       //!< refernce to a grid object.

  //! Enumeration of available output format
  enum class
  output_format : idx_t {
    csv = 0,           //!< comma separated ascii file with headers,
                       //! can be read by common spreadsheet apps
                       //! or by Paraview or Octave.
      
    octave_ascii = 1,  //!< GNU Octave ascii data format, can be
                       //! loaded via the `load` command in GNU Octave.

    json = 2           //! JSON ascii data format, can be reas back in
                       //! via the ctor, useful for restart data
  };

  //! @brief The default generator function used to set up
  //! x-coordinates of particle positions if none is
  //! is specified.

  //! Generates a uniform random distribution.
  double
  default_x_generator ();

  //! @brief The default generator function used to set up
  //! y-coordinates of particle positions if none is
  //! is specified.

  //! Generates a uniform random distribution.
  double
  default_y_generator ();

  //! @brief Template for export function.

  //! If a format is
  //! not specified, just outputs an error message.
  template<output_format fmt>
  void
  print (std::ostream & os) const {
    os << "output format not implementd" << std::endl;
  }
  
  //! @brief Simplest form of constructor.
  
  //! Particle positions are not assigned, they must be set manually later.
  //! @param n number of particles
  //! @param grid_ quadgrid_t object, sizes need to have been already set up.
  particles_t (idx_t n, const quadgrid_t<std::vector<double>>& grid_)
    : num_particles(n), grid(grid_) { } 

  //! @brief Ctor to import data from json.
  
  //! Grid data may be stored in the same `json` object but must be read
  //! separately before invoking this constructor.
  particles_t (const nlohmann::json &j,
	       const quadgrid_t<std::vector<double>>& grid_)
    :  grid(grid_)
  {
    j["dprops"].get_to<dprops_t> (dprops);
    j["iprops"].get_to<iprops_t> (iprops);
    j["x"].get_to<std::vector<double>> (x);
    j["y"].get_to<std::vector<double>> (y);
    j["num_particles"].get_to<idx_t> (num_particles);
  }
  
  //! @brief Constructor with default position generators.
  
  //! Distributes particles randomly over the
  //! grid.
  //! @param n number of particles
  //! @param grid_ quadgrid_t object, sizes need to have been already set up.
  //! @param ipropnames keys for entries in the particles_t::iprops map.
  //! @param dpropnames keys for entries in the particles_t::dprops map.
  particles_t (idx_t n, const std::vector<std::string>& ipropnames,
	       const std::vector<std::string>& dpropnames,
	       const quadgrid_t<std::vector<double>>& grid_);

  //! @brief Constructor with custom position vectors.
  
  //! Distributes particles based on the given vectors
  //! `xv` and `yv`. The vectors are copied, and left unchangend
  //! they can be deleted to reclaim memory if needed.
  //! @param n number of particles
  //! @param grid_ quadgrid_t object, sizes need to have been already set up.
  //! @param ipropnames keys for entries in the particles_t::iprops map.
  //! @param dpropnames keys for entries in the particles_t::dprops map.
  //! @param xv the x-coordinates of new particles.
  //! @param yv the y-coordinates of new particles.
  particles_t (idx_t n, const std::vector<std::string>& ipropnames,
	       const std::vector<std::string>& dpropnames,
	       const quadgrid_t<std::vector<double>>& grid_,
	       const std::vector<double> & xv,
	       const std::vector<double> & yv);

  //! @brief Constructor with custom position generators.
  
  //! Distributes particles based on the given generator functions
  //! `xgen` and `ygen`. Each call to these functions should return
  //! the x- and y-coordinate of a new particle.
  //! @param n number of particles
  //! @param grid_ quadgrid_t object, sizes need to have been already set up.
  //! @param ipropnames keys for entries in the particles_t::iprops map.
  //! @param dpropnames keys for entries in the particles_t::dprops map.
  //! @param xgen generator function for the x-coordinates of new particles.
  //! @param ygen generator function for the y-coordinates of new particles.
  particles_t (idx_t n, const std::vector<std::string>& ipropnames,
	       const std::vector<std::string>& dpropnames,
	       const quadgrid_t<std::vector<double>>& grid_,
	       std::function<double ()> xgen,
	       std::function<double ()> ygen);

  //! @brief Initialize particle properties.
  
  //! Allocates vectors to store particle properties, this is
  //! invoked automatically if the CTOR is invoked specifying
  //! property names, must be invoked manually otherwise.
  //! @param ipropnames keys for entries in the particles_t::iprops map.
  //! @param dpropnames keys for entries in the particles_t::dprops map.
  void
  init_props (const std::vector<std::string>& ipropnames,
	      const std::vector<std::string>& dpropnames);

  //! @brief Erase particcles based on coordinates.

  //! Given a function to decide whether a particle
  //! lies inside a region or not, remove all particles
  //! for which the function returns true, and also 
  //! erase corresponding entries in dprops and iprops.
  void
  remove_in_region (std::function<bool (double, double)> fun) {
    std::vector<idx_t> vin{};
    for (idx_t i = 0; i < num_particles; ++i) {
      if (fun (x[i], y[i])) {
	vin.push_back (i);
      }
    }

    for (auto &dprop : dprops) {
      for (auto id = vin.rbegin (); id != vin.rend (); ++id) {
	dprop.second.erase (dprop.second.begin () + (*id));
      }
    }

    for (auto &iprop : iprops) {
      for (auto id = vin.rbegin (); id != vin.rend (); ++id) {
	iprop.second.erase (iprop.second.begin () + (*id));
      }
    }

    for (auto id = vin.rbegin (); id != vin.rend (); ++id) {
      x.erase (x.begin () + (*id));
      y.erase (y.begin () + (*id));
      --num_particles;
    }
  };

  //! @brief Build grid/particles connectivity.
  
  //! Builds/updates the `grd_to_ptcl` map. Must be used whenever
  //! particles cross cell boundaries.
  void
  init_particle_mesh ();

  //! @brief Initialize particle positions with generator functions.
  
  //! Invoked automatically if the generators are passed to the CTOR,
  //! must be invoked manually otherwise.
  void
  init_particle_positions (std::function<double ()> xgentr,
			   std::function<double ()> ygentr);

  //! @brief Construct a mass matrix.

  //! Must be invoked manually before invoking any of the transfer
  //! methods with flag `use_mass` set to `true`
  void
  build_mass ();

  //! @brief shortcut for `dprops.at (name) [ii]`
  double &
  dp (const std::string & name, idx_t ii) {
    return dprops.at (name) [ii];
  }

  //! @brief shortcut for `dprops.at (name) [ii]`
  const double &
  dp (const std::string & name, idx_t ii) const {
    return dprops.at (name) [ii];
  }
  
  //! @brief shortcut for `iprops.at (name) [ii]`
  idx_t &
  ip (const std::string & name, idx_t ii) {
    return iprops.at (name) [ii];
  }

  //! @brief shortcut for `iprops.at (name) [ii]`
  const idx_t &
  ip (const std::string & name, idx_t ii) const {
    return iprops.at (name) [ii];
  }

  static
  const std::string &
  getkey(std::map<std::string, std::vector<double>> const &varnames,
	 std::size_t ivar)  {
    return std::next (varnames.begin (), ivar)->first;
  };

  static
  const std::string &
  getkey(std::vector<std::string> const &varnames,
	 std::size_t ivar)  {
    return varnames[ivar];
  };

  static
  const char*
  getkey(std::initializer_list<const char *> const &varnames,
	 std::size_t ivar)  {
    return *(std::next (varnames.begin (), ivar));
  };

  //! @brief Map particle variables to the grid.

  //! Assume all fields of `vars` are to be mapped,
  //! and use the same field names for particle and
  //! grid variables.
  void
  p2g (std::map<std::string, std::vector<double>> & vars,
       bool apply_mass = false) const {
    p2g (vars, vars, vars, apply_mass);
  }

  //! @brief Map particle variables to the grid.

  //! Choose which quantities need to be mapped according
  //! to the strings in `gvarnames`,
  //! and use the same field names for particle and
  //! grid variables.
  template<typename GT, typename PT>
  void
  p2g (std::map<std::string, std::vector<double>> & vars,
       PT const & pvarnames,
       GT const & gvarnames,
       bool apply_mass = false,
       assignment_t OP = ASSIGNMENT_OPS::PLUS_EQ) const;

  template<typename str>
  void
  p2g (std::map<std::string, std::vector<double>> & vars,
       std::initializer_list<str> const & pvarnames,
       std::initializer_list<str> const & gvarnames,
       bool apply_mass = false,
       assignment_t OP = ASSIGNMENT_OPS::PLUS_EQ) const;
  
  template<typename GT, typename PT>
  void
  p2gd (std::map<std::string, std::vector<double>> & vars,
	PT const & pxvarnames,
	PT const & pyvarnames,
	std::string const &area,
	GT const & gvarnames,
	bool apply_mass = false,
	assignment_t OP = ASSIGNMENT_OPS::PLUS_EQ) const;

  template<typename str>
  void
  p2gd (std::map<std::string, std::vector<double>> & vars,
	std::initializer_list<str> const & pxvarnames,
	std::initializer_list<str> const & pyvarnames,
	std::string const & area,
	std::initializer_list<str> const & gvarnames,
	bool apply_mass = false,
	assignment_t OP = ASSIGNMENT_OPS::PLUS_EQ) const;

  void
  g2p (const std::map<std::string, std::vector<double>>& vars,
       bool apply_mass = false, assignment_t OP = ASSIGNMENT_OPS::PLUS_EQ) {
    g2p (vars, vars, vars, apply_mass, OP);
  }

  template<typename str>
  void
  g2p (const std::map<std::string, std::vector<double>>& vars,
       std::initializer_list<str> const & gvarnames,
       std::initializer_list<str> const & pvarnames,
       bool apply_mass = false,
       assignment_t OP = ASSIGNMENT_OPS::PLUS_EQ);

  template<typename GT, typename PT>
  void
  g2p (const std::map<std::string, std::vector<double>>& vars,
       GT const & gvarnames,
       PT const & pvarnames,
       bool apply_mass = false,
       assignment_t OP = ASSIGNMENT_OPS::PLUS_EQ);

  template<typename GT, typename PT>
  void
  g2pd (const std::map<std::string, std::vector<double>>& vars,
	GT const & gvarnames,
	PT const & pxvarnames,
	PT const & pyvarnames,
	bool apply_mass = false,
	assignment_t OP = ASSIGNMENT_OPS::PLUS_EQ);

  template<typename str>
  void
  g2pd (const std::map<std::string, std::vector<double>>& vars,
	std::initializer_list<str> const & gvarnames,
	std::initializer_list<str> const &pxvarnames,
	std::initializer_list<str> const & pyvarnames,
	bool apply_mass = false,
	assignment_t OP = ASSIGNMENT_OPS::PLUS_EQ);

};

//! @brief Adaptor to allow implicit conversion from
//! `particles_t` to `json`.
void
to_json (nlohmann::json &j, const particles_t &p);


#include "particles_imp.h"

#endif /* PARTICLES_H */
