#ifndef PARTICLES_H
#define PARTICLES_H

#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <quadgrid_config.h>
#include <json.hpp>
#include <map>
#include <quadgrid_cpp.h>
#include <string>

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
  using idx_t = quadgrid_t<vector_t<real_t>>::idx_t;

  idx_t num_particles;    //!< number of particles.
  vector_t<real_t> x;  //!< x coordinate of particle positions.
  vector_t<real_t> y;  //!< y coordinate of particle positions.

  //! integer type quantities associated with the particles.
  std::map<std::string, vector_t<idx_t>> iprops;

  //! `double` type quantities associated with the particles.
  std::map<std::string, vector_t<real_t>> dprops;

  vector_t<real_t> M; //!< Mass matrix to be used for transfers if required.
  std::map<idx_t, vector_t<idx_t>> grd_to_ptcl;   //!< grid->particles connectivity.
  vector_t<idx_t> ptcl_to_grd;                    //!< particles->grid connectivity.
  vector_t<idx_t> ptcl_grd_color;                 //!< color of particle's cell.
  const quadgrid_t<vector_t<real_t>>& grid;       //!< refernce to a grid object.

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

  //! Enumeration of cell colors
  enum class
  cell_color : idx_t {
    red   = 0,
    green = 1,
    blue  = 2,
    black = 3,
    all_colors = 4
  };

  
  //! @brief The default generator function used to set up
  //! x-coordinates of particle positions if none is
  //! is specified.

  //! Generates a uniform random distribution.
  real_t
  default_x_generator ();

  //! @brief The default generator function used to set up
  //! y-coordinates of particle positions if none is
  //! is specified.

  //! Generates a uniform random distribution.
  real_t
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
  particles_t (idx_t n, const quadgrid_t<vector_t<real_t>>& grid_)
    : num_particles(n), grid(grid_) { }

  //! @brief Ctor to import data from json.

  //! Grid data may be stored in the same `json` object but must be read
  //! separately before invoking this constructor.
  particles_t (const nlohmann::json &j,
               const quadgrid_t<vector_t<real_t>>& grid_)
    :  grid(grid_)
  {
    j["dprops"].get_to<std::map<std::string, vector_t<real_t>>> (dprops);
    j["iprops"].get_to<std::map<std::string, vector_t<int>>> (iprops);
    j["x"].get_to<vector_t<real_t>> (x);
    j["y"].get_to<vector_t<real_t>> (y);
    j["num_particles"].get_to<idx_t> (num_particles);
  }

  //! @brief Constructor with default position generators.

  //! Distributes particles randomly over the
  //! grid.
  //! @param n number of particles
  //! @param grid_ quadgrid_t object, sizes need to have been already set up.
  //! @param ipropnames keys for entries in the particles_t::iprops map.
  //! @param dpropnames keys for entries in the particles_t::dprops map.
  particles_t (idx_t n, const vector_t<std::string>& ipropnames,
               const vector_t<std::string>& dpropnames,
               const quadgrid_t<vector_t<real_t>>& grid_);

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
  particles_t (idx_t n, const vector_t<std::string>& ipropnames,
               const vector_t<std::string>& dpropnames,
               const quadgrid_t<vector_t<real_t>>& grid_,
               const vector_t<real_t> & xv,
               const vector_t<real_t> & yv);

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
  particles_t (idx_t n, const vector_t<std::string>& ipropnames,
               const vector_t<std::string>& dpropnames,
               const quadgrid_t<vector_t<real_t>>& grid_,
               std::function<real_t ()> xgen,
               std::function<real_t ()> ygen);

  //! @brief Initialize particle properties.

  //! Allocates vectors to store particle properties, this is
  //! invoked automatically if the CTOR is invoked specifying
  //! property names, must be invoked manually otherwise.
  //! @param ipropnames keys for entries in the particles_t::iprops map.
  //! @param dpropnames keys for entries in the particles_t::dprops map.
  void
  init_props (const vector_t<std::string>& ipropnames,
              const vector_t<std::string>& dpropnames);

  //! @brief Erase particcles based on coordinates.

  //! Given a function to decide whether a particle
  //! lies inside a region or not, remove all particles
  //! for which the function returns true, and also
  //! erase corresponding entries in dprops and iprops.
  void
  remove_in_region (std::function<bool (real_t, real_t)> fun) {
    vector_t<idx_t> vin{};
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

  //! Builds/updates the `grd_to_ptcl` and `ptcl_to_grd` maps.
  //! Must be used whenever particles cross cell boundaries.
  void
  init_particle_mesh ();

  //! Updates the`ptcl_to_grd` map only, without changing
  //  `grd_to_ptcl`, if not needed.
  void
  update_ptcl_to_grd ();
  
  //! @brief Mark particles by cell color

  void
  mark_by_cell_color ();

  //! @brief Reorder coordinates an properties according to the ordering vvector.
  void
  reorder (vector_t<idx_t> &);

  //! @brief Initialize particle positions with generator functions.

  //! Invoked automatically if the generators are passed to the CTOR,
  //! must be invoked manually otherwise.
  void
  init_particle_positions (std::function<real_t ()> xgentr,
                           std::function<real_t ()> ygentr);

  //! @brief Construct a mass matrix.

  //! Must be invoked manually before invoking any of the transfer
  //! methods with flag `use_mass` set to `true`
  void
  build_mass ();

  //! @brief shortcut for `dprops.at (name) [ii]`
  real_t &
  dp (const std::string & name, idx_t ii) {
    return dprops.at (name) [ii];
  }

  //! @brief shortcut for `dprops.at (name) [ii]`
  const real_t &
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
  getkey(std::map<std::string, vector_t<real_t>> const &varnames,
         std::size_t ivar)  {
    return std::next (varnames.begin (), ivar)->first;
  };

  static
  const std::string &
  getkey(vector_t<std::string> const &varnames,
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
  p2g (std::map<std::string, vector_t<real_t>> & vars,
       bool apply_mass = false)  {
    p2g (vars, vars, vars, apply_mass);
  }

  //! @brief Map particle variables to the grid.

  //! Choose which quantities need to be mapped according
  //! to the strings in `gvarnames`,
  //! and use the same field names for particle and
  //! grid variables.
  template<typename GT, typename PT>
  void
  p2g (std::map<std::string, vector_t<real_t>> & vars,
       PT const & pvarnames,
       GT const & gvarnames,
       bool apply_mass = false) ;

  template<typename str>
  void
  p2g (std::map<std::string, vector_t<real_t>> & vars,
       std::initializer_list<str> const & pvarnames,
       std::initializer_list<str> const & gvarnames,
       bool apply_mass = false) ;

  template<typename GT, typename PT>
  void
  p2gd (std::map<std::string, vector_t<real_t>> & vars,
        PT const & pxvarnames,
        PT const & pyvarnames,
        std::string const &area,
        GT const & gvarnames,
        bool apply_mass = false);

  template<typename str>
  void
  p2gd (std::map<std::string, vector_t<real_t>> & vars,
        std::initializer_list<str> const & pxvarnames,
        std::initializer_list<str> const & pyvarnames,
        std::string const & area,
        std::initializer_list<str> const & gvarnames,
        bool apply_mass = false);

  void
  g2p (const std::map<std::string, vector_t<real_t>>& vars,
       bool apply_mass = false) {
    g2p (vars, vars, vars, apply_mass);
  }

  template<typename str>
  void
  g2p (const std::map<std::string, vector_t<real_t>>& vars,
       std::initializer_list<str> const & gvarnames,
       std::initializer_list<str> const & pvarnames,
       bool apply_mass = false);

  template<typename GT, typename PT>
  void
  g2p (const std::map<std::string, vector_t<real_t>>& vars,
       GT const & gvarnames,
       PT const & pvarnames,
       bool apply_mass = false);

  template<typename GT, typename PT>
  void
  g2pd (const std::map<std::string, vector_t<real_t>>& vars,
        GT const & gvarnames,
        PT const & pxvarnames,
        PT const & pyvarnames,
        bool apply_mass = false);

  template<typename str>
  void
  g2pd (const std::map<std::string, vector_t<real_t>>& vars,
        std::initializer_list<str> const & gvarnames,
        std::initializer_list<str> const &pxvarnames,
        std::initializer_list<str> const & pyvarnames,
        bool apply_mass = false);

};

//! @brief Adaptor to allow implicit conversion from
//! `particles_t` to `json`.
void
to_json (nlohmann::json &j, const particles_t &p);

//! @brief Template class for the update of 
//! `ptcl_to_grd` mapping.
template<typename P2G_t, typename COORD_t>
class
ptcl_to_grd_update_t {

  using idx_t=particles_t::idx_t;
  P2G_t ptcl_to_grd;
  const COORD_t x;
  const COORD_t y;
  const real_t hx;
  const real_t hy;
  const idx_t nrows;
  
public :
  ptcl_to_grd_update_t (P2G_t ptcl_to_grd_,
			const COORD_t x_, const COORD_t y_,
			real_t hx_, real_t hy_, const idx_t nrows_)
    : ptcl_to_grd(ptcl_to_grd_), x(x_), y(y_), hx(hx_), hy(hy_), nrows(nrows_) { }

  void operator() (particles_t::idx_t ii) {
    ptcl_to_grd[ii] = quadgrid_t<COORD_t>::sub2gind (static_cast<idx_t> (std::floor (y[ii] / hy)),
						     static_cast<idx_t> (std::floor (x[ii] / hx)),
						     nrows);
  }
};


#include "particles_imp.h"

#endif /* PARTICLES_H */
