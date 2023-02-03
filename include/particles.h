#ifndef PARTICLES_H
#define PARTICLES_H

#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <quadgrid_cpp.h>
#include <string>


//! \brief Class to represent particles embedded in a grid.

//! Offers methods for transfer of quantities from particles
//! to grid or vice-versa (`p2g`, `g2p, p2gd`, `g2pd`).
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
  using idx_t = quadgrid_t<std::vector<double>>::idx_t;
  
  idx_t num_particles;    //!< number of particles.
  std::vector<double> x;  //!< x coordinate of particle positions.
  std::vector<double> y;  //!< y coordinate of particle positions.

  std::map<std::string, std::vector<idx_t>> iprops;   //!< integer type quantities associated with the particles.
  std::map<std::string, std::vector<double>> dprops;  //!< `double` type quantities associated with the particles.

  std::vector<double> M;                             //!< Mass matrix to be used for transfers if required.
  std::map<idx_t, std::vector<idx_t>> grd_to_ptcl;   //!< grid/particles connectivity.
  const quadgrid_t<std::vector<double>>& grid;       //!< refernce to a grid object.

  //! Enumeration of available output format
  enum class
  output_format : idx_t {
    csv = 0,           //!< comma separated ascii file with headers, can be read by common spreadsheet apps or by Paraview or Octave.
    octave_ascii = 1   //!< GNU Octave ascii data format, can be loaded via the `load` command in GNU Octave.
  };

  double
  default_x_generator ();

  double
  default_y_generator ();

  template<output_format fmt = output_format::csv>
  void
  print (std::ostream & os) const {
    os << "output format not implementd" << std::endl;
  }

  particles_t (idx_t n, const quadgrid_t<std::vector<double>>& grid_)
    : num_particles(n), grid(grid_) { }


  particles_t (idx_t n, const std::vector<std::string>& ipropnames,
	       const std::vector<std::string>& dpropnames,
	       const quadgrid_t<std::vector<double>>& grid_);

  particles_t (idx_t n, const std::vector<std::string>& ipropnames,
	       const std::vector<std::string>& dpropnames,
	       const quadgrid_t<std::vector<double>>& grid_,
	       const std::vector<double> & xgen,
	       const std::vector<double> & ygen);

  particles_t (idx_t n, const std::vector<std::string>& ipropnames,
	       const std::vector<std::string>& dpropnames,
	       const quadgrid_t<std::vector<double>>& grid_,
	       std::function<double ()> xgen,
	       std::function<double ()> ygen);

  void
  init_props (const std::vector<std::string>& ipropnames,
	      const std::vector<std::string>& dpropnames);

  void
  init_particle_mesh ();
  
  void
  init_particle_positions (std::function<double ()> xgentr,
			   std::function<double ()> ygentr);
  
  void
  build_mass ();

  const std::string &
  getkey(std::map<std::string, std::vector<double>> const &varnames,
	 int ivar) const {
    return std::next (varnames.begin (), ivar)->first;
  };


  const std::string &
  getkey(std::vector<std::string> const &varnames,
	 std::size_t ivar) const {
    return varnames[ivar];
  };

  void
  p2g (std::map<std::string, std::vector<double>> & vars,
       bool apply_mass = false) const {
    p2g (vars, vars, vars, apply_mass);
  }

  template<typename GT, typename PT>
  void
  p2g (std::map<std::string, std::vector<double>> & vars,
       PT const & pvarnames,
       GT const & gvarnames,
       bool apply_mass = false) const;

  template<typename GT, typename PT>
  void
  p2gd (std::map<std::string, std::vector<double>> & vars,
	PT const & pxvarnames,
	PT const & pyvarnames,
	std::string const &area,
	GT const & gvarnames,
	bool apply_mass = false) const;

  void
  g2p (const std::map<std::string, std::vector<double>>& vars,
       bool apply_mass = false) {
    g2p (vars, vars, vars, apply_mass);
  }

  template<typename GT, typename PT>
  void
  g2p (const std::map<std::string, std::vector<double>>& vars,
       GT const & gvarnames,
       PT const & pvarnames,
       bool apply_mass = false);

  template<typename GT, typename PT>
  void
  g2pd (const std::map<std::string, std::vector<double>>& vars,
	GT const & gvarnames,
	PT const & pxvarnames,
	PT const & pyvarnames,
	bool apply_mass = false);

};

#include "particles_imp.h"

#endif /* PARTICLES_H */
