#include <algorithm>
#include <iomanip>
#include <iostream>
#include <random>

#include <particles.h>

namespace ASSIGNMENT_OPS {
  assignment_t EQ = [] (double& TO, const double& FROM) -> double& { return TO = FROM; };
  assignment_t PLUS_EQ = [] (double& TO, const double& FROM) -> double& { return TO += FROM; };
  assignment_t TIMES_EQ = [] (double& TO, const double& FROM) -> double& { return TO *= FROM; };
}


double
particles_t::default_x_generator () {
  static std::random_device rd;
  static std::mt19937 gen (rd ());
  static std::uniform_real_distribution<> dis (0.0, 1.0);
  return dis (gen) * grid.num_cols () * grid.hx ();
}

double
particles_t::default_y_generator () {
  static std::random_device rd;
  static std::mt19937 gen (rd ());
  static std::uniform_real_distribution<> dis (0.0, 1.0);
  return dis (gen) * grid.num_rows () * grid.hy ();
}

particles_t::particles_t
(
 idx_t n, const std::vector<std::string>& ipropnames,
 const std::vector<std::string>& dpropnames,
 const quadgrid_t<std::vector<double>>& grid_
 ) :  particles_t (n, grid_) {

  init_props (ipropnames, dpropnames);

  init_particle_positions
    (
     [this] { return this->default_x_generator (); },
     [this] { return this->default_y_generator (); }
     );

  init_particle_mesh ();
}


particles_t::particles_t
(
 idx_t n, const std::vector<std::string>& ipropnames,
 const std::vector<std::string>& dpropnames,
 const quadgrid_t<std::vector<double>>& grid_,
 const std::vector<double> & xgen,
 const std::vector<double> & ygen
 ) : particles_t (n, grid_) {

  x = xgen;
  y = ygen;

  init_props (ipropnames, dpropnames);

  init_particle_mesh ();
}

particles_t::particles_t
(
 idx_t n, const std::vector<std::string>& ipropnames,
 const std::vector<std::string>& dpropnames,
 const quadgrid_t<std::vector<double>>& grid_,
 std::function<double ()> xgen,
 std::function<double ()> ygen
 ) : particles_t (n, grid_) {

  init_props (ipropnames, dpropnames);

  init_particle_positions (xgen, ygen);

  init_particle_mesh ();
}


void
particles_t::init_props
(
 const std::vector<std::string>& ipropnames,
 const std::vector<std::string>& dpropnames
 ) {

  for (idx_t ii = 0; ii < ipropnames.size (); ++ii) {
    iprops[ipropnames[ii]].assign (num_particles, 0);
  }

  for (idx_t ii = 0; ii < dpropnames.size (); ++ii) {
    dprops[dpropnames[ii]].assign (num_particles, 0.0);
  }
}


void
particles_t::init_particle_mesh () {

  for (auto & igrd : grd_to_ptcl)
    std::vector<idx_t>{}. swap (igrd.second);

  for (auto ii = 0; ii < x.size (); ++ii) {
    idx_t c = static_cast<idx_t> (std::floor (x[ii] / grid.hx ()));
    idx_t r = static_cast<idx_t> (std::floor (y[ii] / grid.hy ()));

    grd_to_ptcl[grid.sub2gind (r, c)].push_back (ii);
  }
}


void
particles_t::init_particle_positions
(
 std::function<double ()> xgentr,
 std::function<double ()> ygentr
 ) {
  x.resize (num_particles);
  y.resize (num_particles);

  std::generate (x.begin (), x.end (), xgentr);
  std::generate (y.begin (), y.end (), ygentr);
}



void
particles_t::build_mass () {
  M.assign (grid.num_global_nodes (), 0.0);
  for (auto icell = grid.begin_cell_sweep ();
       icell != grid.end_cell_sweep (); ++icell) {
    for (auto inode = 0;
         inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell;
         ++inode) {
      M[icell->gt (inode)] += (grid.hx () / 2.) * (grid.hy () / 2.);
    }
  }
}

template<>
void
particles_t::print<particles_t::output_format::json> (std::ostream & os) const {
  nlohmann::json j = *this;
  os << j;
}

template<>
void
particles_t::print<particles_t::output_format::csv> (std::ostream & os) const {

  os << "\"x\", " << "\"y\"";

  for (auto const & ii : dprops)
    os << ", \"" << ii.first << "\"";

  for (auto const & ii : iprops)
    os << ", \"" << ii.first << "\"";

  os << std::endl;


  for (idx_t jj = 0; jj < x.size (); ++jj) {

    os << x[jj] << ", " << y[jj];

    for (auto const & ii : dprops)
      os << ", " << std::setprecision (16) << ii.second[jj];

    for (auto const & ii : iprops)
      os << ", " << ii.second[jj];

    os << std::endl;
  }
}

template<>
void
particles_t::print<particles_t::output_format::octave_ascii>
(std::ostream & os) const {

  os << "# name: x" << std::endl
     << "# type: matrix" << std::endl
     << "# rows: 1" << std::endl
     << "# columns: " << x.size () << std::endl;
  for (auto const & kk : x) {
    os  << std::setprecision(16) << kk << " ";
  }
  os << std::endl;

  os << "# name: y" << std::endl
     << "# type: matrix" << std::endl
     << "# rows: 1" << std::endl
     << "# columns: " << y.size () << std::endl;
  for (auto const & kk : y) {
    os  << std::setprecision(16) << kk << " ";
  }
  os << std::endl;

  os << "# name: dprops" << std::endl
     << "# type: scalar struct" << std::endl
     << "# ndims: 2" << std::endl
     << "1 1" << std::endl
     << "# length: " << dprops.size () << std::endl;


  for (auto const & ii : dprops) {
    os << "# name: " << ii.first << std::endl
       << "# type: matrix" << std::endl
       << "# rows: 1" << std::endl
       << "# columns: " << ii.second.size () << std::endl;
    for (auto const & kk : ii.second) {
      os << std::setprecision(16) << kk << " ";
    }
    os << std::endl;
  }
  os << std::endl;

  os << "# name: iprops" << std::endl
     << "# type: scalar struct" << std::endl
     << "# ndims: 2" << std::endl
     << "1 1" << std::endl
     << "# length: " << iprops.size () << std::endl;

  for (auto const & ii : iprops) {
    os << "# name: " << ii.first << std::endl
       << "# type: int64 matrix" << std::endl
       << "# ndims: 2" << std::endl
       << "1 " << ii.second.size () << std::endl;
    for (auto const & kk : ii.second) {
      os << kk << " ";
    }
    os << std::endl;
  }
  os << std::endl;
}

void
to_json (nlohmann::json &j, const particles_t &p) {
  j = nlohmann::json{
    {"num_particles", p.num_particles},
    {"x", p.x},
    {"y", p.y},
    {"dprops", p.dprops},
    {"iprops", p.iprops},
    {"grid_properties",
     {{"nx", p.grid.num_cols ()},
      {"ny", p.grid.num_rows ()},
      {"hx", p.grid.hx ()},
      {"hy", p.grid.hy ()}}
    }
  };
}

void
particles_t::reorder (std::vector<idx_t> &ordering) {

  for (idx_t ii = 0; ii < num_particles - 1; ++ii) {
    if (ii != ordering[ii]) {

      for (auto &dprop : dprops) {
        auto &col = dprop.second;
        std::swap (col[ii], col[ordering[ii]]);
      }

      for (auto &iprop : iprops) {
        auto &col = iprop.second;
        std::swap (col[ii], col[ordering[ii]]);
      }

      std::swap (x[ii], x[ordering[ii]]);
      std::swap (y[ii], y[ordering[ii]]);

      for (int jj = ii; jj < ordering.size (); ++jj) {
        if (ordering[jj] == ii) {
          ordering[jj] = ordering[ii];
          ordering[ii] = ii;
          break;
        }
      }

    }
  }

}
