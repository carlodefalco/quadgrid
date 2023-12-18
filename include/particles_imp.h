#ifndef PARTICLES_IMP_H
#define PARTICLES_IMP_H

template<typename str>
void
particles_t::p2g
(std::map<std::string, std::vector<double>> & vars,
 std::initializer_list<str> const & pvarnames,
 std::initializer_list<str> const & gvarnames,
 bool apply_mass) const {
  using strlist = std::initializer_list<str> const &;
  p2g<strlist, strlist>
    (vars, pvarnames, gvarnames, apply_mass);
}

template<typename GT, typename PT>
void
particles_t::p2g
(std::map<std::string, std::vector<double>> & vars,
 PT const & pvarnames,
 GT const & gvarnames,
 bool apply_mass) const {

  using idx_t = quadgrid_t<std::vector<double>>::idx_t;
  double N = 0.0, xx = 0.0, yy = 0.0;
  idx_t idx = 0;


  for (std::size_t ivar = 0; ivar < std::size(gvarnames); ++ivar) {
    auto & gvar = vars[getkey(gvarnames, ivar)];
    auto const & dprop = dprops.at (getkey(pvarnames, ivar));

    // for (auto icell = grid.begin_cell_sweep ();
    //   icell != grid.end_cell_sweep (); ++icell) {
    //   if (grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
    //  for (idx_t ii = 0;
    //       ii < grd_to_ptcl.at (icell->get_global_cell_idx ()).size ();
    //       ++ii) {
    //    idx = grd_to_ptcl.at (icell->get_global_cell_idx ())[ii];
    //    xx = x[idx];
    //    yy = y[idx];
    //    for (idx_t inode = 0; inode < 4; ++inode) {
    //      N = icell->shp(xx, yy, inode);
    //      gvar[icell->gt(inode)] +=
    //        N * dprop[idx];
    //    }
    //  }
    // }


    for (idx_t ip = 0; ip <= this->num_particles; ++ip) {
      xx = x[ip];
      yy = y[ip];
      auto icell = grid[ptcl_to_grd[ip]];
      for (idx_t inode = 0; inode < 4; ++inode) {
        N = icell.shp(xx, yy, inode);

        gvar[icell.gt(inode)] +=
          N * dprop[ip];
      }
    }


  }

  if (apply_mass)
    for (std::size_t ivar = 0; ivar < std::size (gvarnames); ++ivar)
      for (idx_t ii = 0; ii < M.size (); ++ii) {
        vars[getkey(gvarnames, ivar)][ii]  /= M[ii];
      }
}


template<typename str>
void
particles_t::p2gd
(std::map<std::string, std::vector<double>> & vars,
 std::initializer_list<str> const & pxvarnames,
 std::initializer_list<str> const & pyvarnames,
 std::string const &area,
 std::initializer_list<str> const & gvarnames,
 bool apply_mass) const {
  using strlist = std::initializer_list<str> const &;
  p2gd<strlist, strlist>
    (vars, pxvarnames, pyvarnames, area, gvarnames, apply_mass);
}


template<typename GT, typename PT>
void
particles_t::p2gd
(std::map<std::string, std::vector<double>> & vars,
 PT const & pxvarnames,
 PT const & pyvarnames,
 std::string const &area,
 GT const & gvarnames,
 bool apply_mass) const {

  using idx_t = quadgrid_t<std::vector<double>>::idx_t;
  double xx = 0.0, yy = 0.0, Nx = 0.0, Ny = 0.0;
  //idx_t idx = 0;

  for (std::size_t ivar = 0; ivar < std::size (gvarnames); ++ivar) {
    auto & gvar = vars[getkey(gvarnames, ivar)];
    auto const & dpropx = dprops.at (getkey(pxvarnames, ivar));
    auto const & dpropy = dprops.at (getkey(pyvarnames, ivar));
    auto const & dproparea = dprops.at (area);


    // for (auto icell = grid.begin_cell_sweep();
    //   icell != grid.end_cell_sweep(); ++icell) {
    //   if (grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
    //          for (idx_t ii = 0;
    //               ii < grd_to_ptcl.at (icell->get_global_cell_idx ()).size ();
    //               ++ii) {
    //            idx = grd_to_ptcl.at (icell->get_global_cell_idx())[ii];
    //            xx = x[idx];
    //            yy = y[idx];
    //            for (idx_t inode=0; inode<4; ++inode) {
    //              Nx = icell->shg (xx, yy, 0, inode);
    //              Ny = icell->shg (xx, yy, 1, inode);
    //              gvar[icell->gt(inode)] +=
    //                (Nx * dpropx[idx] + Ny * dpropy[idx]) * dproparea[idx];
    //            }
    //          }
    // }

    for (idx_t ip = 0; ip <= this->num_particles; ++ip) {
      xx = x[ip];
      yy = y[ip];
      auto icell = grid[ptcl_to_grd[ip]];
      for (idx_t inode = 0; inode < 4; ++inode) {
        Nx = icell.shg (xx, yy, 0, inode);
        Ny = icell.shg (xx, yy, 1, inode);
        gvar[icell.gt(inode)] += (Nx * dpropx[ip] + Ny * dpropy[ip]) * dproparea[ip];
      }
    }

  }

  if (apply_mass)
    for (std::size_t ivar = 0; ivar < std::size (gvarnames); ++ivar)
      for (idx_t ii = 0; ii < M.size (); ++ii) {
        vars[getkey(gvarnames, ivar)][ii]  /= M[ii];
      }

}

template<typename str>
void
particles_t::g2p
(const std::map<std::string, std::vector<double>> & vars,
 std::initializer_list<str> const & gvarnames,
 std::initializer_list<str> const & pvarnames,
 bool apply_mass) {
  using strlist = std::initializer_list<str> const &;
  g2p<strlist, strlist> (vars, gvarnames,
                         pvarnames, apply_mass);
}


template<typename GT, typename PT>
void
particles_t::g2p
(const std::map<std::string, std::vector<double>>& vars,
 GT const & gvarnames,
 PT const & pvarnames,
 bool apply_mass) {

  using idx_t = particles_t::idx_t;
  double N = 0.0, xx = 0.0, yy = 0.0;
  idx_t idx = 0;

  for (std::size_t ivar = 0; ivar < std::size (gvarnames); ++ivar) {
    auto & dprop = dprops.at (getkey (pvarnames, ivar));
    auto const & gvar = vars.at (getkey (gvarnames, ivar));
    
    // for (auto icell = grid.begin_cell_sweep ();
    //   icell != grid.end_cell_sweep (); ++icell) {
    //   if (grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
    //  for (idx_t ii = 0;
    //       ii < grd_to_ptcl.at (icell->get_global_cell_idx ()).size ();
    //       ++ii) {
    //    idx = grd_to_ptcl.at(icell->get_global_cell_idx ())[ii];
    //    xx = x[idx];
    //    yy = y[idx];
    //    for (idx_t inode = 0; inode < 4; ++inode) {
    //      N = apply_mass ?
    //        icell->shp(xx, yy, inode) * M[icell->gt(inode)] :
    //        icell->shp(xx, yy, inode);
    //      dprop [idx] += N * gvar[icell->gt(inode)];
    //    }
    //  }
    // }

    g2p_helper_t helper (x.begin (), y.begin (), M.cbegin (),
    			 gvar.cbegin (), ptcl_to_grd.cbegin (),
    			 grid.num_rows (), grid.hx (), grid.hy (),
    			 dprop.begin (), apply_mass);
    
    // for (idx_t ip = 0; ip < this->num_particles; ++ip) {
    //   xx = x[ip];
    //   yy = y[ip];
    //   auto icell = grid[ptcl_to_grd[ip]];
    //   for (idx_t inode = 0; inode < 4; ++inode) {
    // 	N = apply_mass ?
    //       icell.shp(xx, yy, inode) * M[icell.gt(inode)] :
    //       icell.shp(xx, yy, inode);
    //     dprop[ip]      += N * gvar[icell.gt(inode)];
    //   }
    // }
       
    range rng (0, this->num_particles);
    std::for_each (rng.begin (), rng.end (), helper);
    
  }
}

template<typename str>
void
particles_t::g2pd
(const std::map<std::string, std::vector<double>>& vars,
 std::initializer_list<str> const & gvarnames,
 std::initializer_list<str> const & pxvarnames,
 std::initializer_list<str> const & pyvarnames,
 bool apply_mass) {
  using strlist = std::initializer_list<str> const &;
  g2pd<strlist, strlist> (vars, gvarnames, pxvarnames,
                          pyvarnames, apply_mass);
}

template<typename GT, typename PT>
void
particles_t::g2pd
(const std::map<std::string, std::vector<double>>& vars,
 GT const & gvarnames,
 PT const & pxvarnames,
 PT const & pyvarnames,
 bool apply_mass) {

  using idx_t = quadgrid_t<std::vector<double>>::idx_t;
  double Nx = 0.0, Ny = 0.0, xx = 0.0, yy = 0.0;
  //idx_t idx = 0;

  for (std::size_t ivar = 0; ivar < std::size (gvarnames); ++ivar) {
    auto const & gvar = vars.at (getkey (gvarnames, ivar));
    auto & dpropx = dprops.at (getkey (pxvarnames, ivar));
    auto & dpropy = dprops.at (getkey (pyvarnames, ivar));

    // for (auto icell = grid.begin_cell_sweep ();
    //   icell != grid.end_cell_sweep (); ++icell) {
    //   if (grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0) {
    //  for (idx_t ii = 0;
    //       ii < grd_to_ptcl.at (icell->get_global_cell_idx ()).size ();
    //       ++ii) {
    //    idx = grd_to_ptcl.at(icell->get_global_cell_idx ())[ii];
    //    xx = x[idx];
    //    yy = y[idx];
    //    for (idx_t inode = 0; inode < 4; ++inode) {
    //      Nx = apply_mass ?
    //        icell->shg(xx, yy, 0, inode) * M[icell->gt(inode)] :
    //        icell->shg(xx, yy, 0, inode);
    //      Ny = apply_mass ?
    //        icell->shg(xx, yy, 1, inode) * M[icell->gt(inode)] :
    //        icell->shg(xx, yy, 1, inode);
    //      dpropx[idx] += Nx * gvar[icell->gt(inode)];
    //      dpropy[idx] += Ny * gvar[icell->gt(inode)];
    //    }
    //  }
    //   }
    // }


    for (idx_t ip = 0; ip <= this->num_particles; ++ip) {
      xx = x[ip];
      yy = y[ip];
      auto icell = grid[ptcl_to_grd[ip]];
      for (idx_t inode = 0; inode < 4; ++inode) {
        Nx = apply_mass ?
          icell.shg(xx, yy, 0, inode) * M[icell.gt(inode)] :
          icell.shg(xx, yy, 0, inode);
        Ny = apply_mass ?
          icell.shg(xx, yy, 1, inode) * M[icell.gt(inode)] :
          icell.shg(xx, yy, 1, inode);
        dpropx[ip] += Nx * gvar[icell.gt(inode)];
        dpropy[ip] += Ny * gvar[icell.gt(inode)];

      }
    }

  }

}

template<>
void
particles_t::print<particles_t::output_format::octave_ascii>
(std::ostream & os) const;

template<>
void
particles_t::print<particles_t::output_format::json>
(std::ostream & os) const;

template<>
void
particles_t::print<particles_t::output_format::csv>
(std::ostream & os) const;

#endif
