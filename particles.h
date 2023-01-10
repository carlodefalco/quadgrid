#ifndef PARTICLES_H
#define PARTICLES_H

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <quadgrid_cpp.h>
#include <random>
#include <string>



struct
particles_t {

  using idx_t = quadgrid_t<std::vector<double>>::idx_t;
  std::vector<double> x;
  std::vector<double> y;

  std::map<std::string, std::vector<idx_t>> iprops;
  std::map<std::string, std::vector<double>> dprops;

  std::vector<double> M;
  std::map<idx_t, std::vector<idx_t>> grd_to_ptcl;
  const quadgrid_t<std::vector<double>>& grid;

  enum class
  output_format : idx_t {
    csv = 0,
    octave_ascii = 1
  };

  double
  default_x_generator () {
    static std::random_device rd;
    static std::mt19937 gen (rd ());
    static std::uniform_real_distribution<> dis (0.0, 1.0);
    return dis (gen) * grid.num_cols () * grid.hx ();
  }
  double
  default_y_generator () {
    static std::random_device rd;
    static std::mt19937 gen (rd ());
    static std::uniform_real_distribution<> dis (0.0, 1.0);
    return dis (gen) * grid.num_rows () * grid.hy (); 
  }
  
  template<output_format fmt = output_format::csv>
  void
  print (std::ostream & os) const {
    os << "output format not implementd" << std::endl;
  };

  template<>
  void
  print<output_format::csv> (std::ostream & os) const {
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
  };

  template<>
  void
  print<output_format::octave_ascii> (std::ostream & os) const {


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
  };

  particles_t (idx_t n, const std::vector<std::string>& ipropnames,
	       const std::vector<std::string>& dpropnames,
	       const quadgrid_t<std::vector<double>>& grid_,
	       std::function<double ()> xgentr = default_x_generator,
	       std::function<double ()> ygentr = default_y_generator)
    : x(n, 0.0), y(n, 0.0), grid(grid_) {

    for (idx_t ii = 0; ii < ipropnames.size (); ++ii) {
      iprops[ipropnames[ii]].assign (n, 0);
    }

    for (idx_t ii = 0; ii < dpropnames.size (); ++ii) {
      dprops[dpropnames[ii]].assign (n, 0.0);
    }

    M = std::vector<double> (grid.num_global_nodes (), 0.0);
    build_mass ();

    init_particle_positions (xgentr, ygentr);

    init_particle_mesh ();
  };

  void
  init_particle_mesh () {

    for (auto & igrd : grd_to_ptcl)
      std::vector<idx_t>{}. swap (igrd.second);

    for (auto ii = 0; ii < x.size (); ++ii) {
      idx_t c = static_cast<idx_t> (std::floor (x[ii] / grid.hx ()));
      idx_t r = static_cast<idx_t> (std::floor (y[ii] / grid.hy ()));

      grd_to_ptcl[grid.sub2gind (r, c)].push_back (ii);
    }
  };

  void
  init_particle_positions (std::function<double ()> xgentr,
			   std::function<double ()> ygentr)
  {
    std::generate (x.begin (), x.end (), xgentr);
    std::generate (y.begin (), y.end (), ygentr);
  };

  void
  build_mass () {
    M.assign (M.size (), 0.0);
    for (auto icell = grid.begin_cell_sweep ();
	 icell != grid.end_cell_sweep (); ++icell) {
      for (auto inode = 0;
	   inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell;
	   ++inode) {
	M[icell->gt (inode)] += (grid.hx () / 2.) * (grid.hy () / 2.);
      }
    }
  };


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
       bool apply_mass = false) const {

    double N = 0.0, xx = 0.0, yy = 0.0;
    idx_t idx = 0;

    for (auto icell = grid.begin_cell_sweep ();
	 icell != grid.end_cell_sweep (); ++icell) {

      if (grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
	for (idx_t ii = 0;
	     ii < grd_to_ptcl.at (icell->get_global_cell_idx ()).size ();
	     ++ii) {
	  idx = grd_to_ptcl.at (icell->get_global_cell_idx ())[ii];
	  xx = x[idx];
	  yy = y[idx];

	  for (idx_t inode = 0; inode < 4; ++inode) {
	    N = icell->shp(xx, yy, inode);
	    for (std::size_t ivar = 0; ivar < gvarnames.size (); ++ivar) {
	      vars[getkey(gvarnames, ivar)][icell->gt(inode)]  +=
		N * dprops.at (getkey(pvarnames, ivar))[idx];
	    }
	  }
	}
    }

    if (apply_mass)
      for (std::size_t ivar = 0; ivar < gvarnames.size (); ++ivar)
	for (idx_t ii = 0; ii < M.size (); ++ii) {
	  vars[getkey(gvarnames, ivar)][ii]  /= M[ii];
	}
  };

  template<typename GT, typename PT>
  void
  p2gd (std::map<std::string, std::vector<double>> & vars,
	PT const & pvarnames,
	GT const & gxvarnames,
	GT const & gyvarnames,
	bool apply_mass = false) const {

    double xx = 0.0, yy = 0.0, Nx = 0.0, Ny = 0.0;
    idx_t idx = 0;

    for (auto icell = grid.begin_cell_sweep();
	 icell != grid.end_cell_sweep(); ++icell) {

      if (grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
	for (idx_t ii = 0;
	     ii < grd_to_ptcl.at (icell->get_global_cell_idx ()).size ();
	     ++ii) {

	  idx = grd_to_ptcl.at (icell->get_global_cell_idx())[ii];
	  xx = x[idx];
	  yy = y[idx];

	  for (idx_t inode=0; inode<4; ++inode) {

	    Nx = icell->shg (xx, yy, 0, inode);
	    Ny = icell->shg (xx, yy, 1, inode);

	    for (std::size_t ivar = 0; ivar < pvarnames.size (); ++ivar) {
	      vars[getkey(gxvarnames, ivar)][icell->gt(inode)]  +=
		Nx * dprops.at (getkey(pvarnames, ivar))[idx];
	      vars[getkey(gyvarnames, ivar)][icell->gt(inode)]  +=
		Ny * dprops.at (getkey(pvarnames, ivar))[idx];
	    }
	  }
	}

    }

    if (apply_mass)
      for (std::size_t ivar = 0; ivar < pvarnames.size (); ++ivar)
	for (idx_t ii = 0; ii < M.size (); ++ii) {
	  vars[getkey(gxvarnames, ivar)][ii]  /= M[ii];
	  vars[getkey(gyvarnames, ivar)][ii]  /= M[ii];
	}

  };

  void
  g2p (const std::map<std::string, std::vector<double>>& vars,
       bool apply_mass) {
    g2p (vars, vars, vars, apply_mass);
  }

  template<typename GT, typename PT>
  void
  g2p (const std::map<std::string, std::vector<double>>& vars,
       GT const & gvarnames,
       PT const & pvarnames,
       bool apply_mass) {

    double N = 0.0, xx = 0.0, yy = 0.0;
    idx_t idx = 0;

    for (auto icell = grid.begin_cell_sweep ();
	 icell != grid.end_cell_sweep (); ++icell) {

      if (grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
	for (idx_t ii = 0;
	     ii < grd_to_ptcl.at (icell->get_global_cell_idx ()).size ();
	     ++ii) {

	  idx = grd_to_ptcl.at(icell->get_global_cell_idx ())[ii];
	  xx = x[idx];
	  yy = y[idx];

	  for (idx_t inode = 0; inode < 4; ++inode) {
	    N = apply_mass ?
	      icell->shp(xx, yy, inode) * M[icell->gt(inode)] :
	      icell->shp(xx, yy, inode);
	    for (std::size_t ivar = 0; ivar < gvarnames.size (); ++ivar)
	      dprops.at (getkey (pvarnames, ivar))[idx] += N * vars.at (getkey (gvarnames, ivar))[icell->gt(inode)];
	  }
	}

    }

  };


};

#endif /* PARTICLES_H */
