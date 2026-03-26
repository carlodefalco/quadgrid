#ifndef PARTICLES_IMP_H
#define PARTICLES_IMP_H
#include "counter.h"

template<typename str>
void
particles_t::p2g
(std::map<std::string, vector_t<real_t>> & vars,
 std::initializer_list<str> const & pvarnames,
 std::initializer_list<str> const & gvarnames,
 bool apply_mass)  {
  using strlist = std::initializer_list<str> const &;
  p2g<strlist, strlist>
    (vars, pvarnames, gvarnames, apply_mass);
}

template<typename GVAR_t, typename PVAR_t, typename P2C_t>
class
p2g_helper_t{

  using idx_t = particles_t::idx_t;
  const PVAR_t x;
  const PVAR_t y;
  const PVAR_t dprop;
  const P2C_t ptcl_to_grd;
  const idx_t nrows;
  const real_t hx;
  const real_t hy;
  GVAR_t M;
  GVAR_t gvar;

  bool apply_mass;

public :

  p2g_helper_t (const PVAR_t x_, const PVAR_t y_, GVAR_t M_, 
		GVAR_t gvar_, const P2C_t ptcl_to_grd_, const idx_t nrows_,
		const real_t hx_, const real_t hy_, const PVAR_t dprop_, bool apply_mass_)
    : x(x_), y(y_), M(M_), gvar(gvar_), 
      ptcl_to_grd(ptcl_to_grd_), nrows(nrows_), hx(hx_), hy(hy_),
      dprop(dprop_), apply_mass(apply_mass_) {};
  
  DEVICE
  void
  operator() (idx_t ip) {
    using qgt = quadgrid_t<GVAR_t>;
    real_t N = 0.0;
    auto xx = x[ip];
    auto yy = y[ip];
    auto r = qgt::gind2row (ptcl_to_grd[ip], nrows);
    auto c = qgt::gind2col (ptcl_to_grd[ip], nrows);
    for (idx_t inode = 0; inode < 4; ++inode) {  
      N = apply_mass ? qgt::shp (xx, yy, inode, c, r, hx, hy) * M[qgt::gt(inode, c, r, nrows)] :
	    qgt::shp (xx, yy, inode, c, r, hx, hy);
      atomicAdd(&(gvar[qgt::gt(inode, c, r, nrows)]), N*dprop[ip]);
    }
  } 


};


template<typename GT, typename PT>
void
particles_t::p2g
(std::map<std::string, vector_t<real_t>> & vars,
 PT const & pvarnames,
 GT const & gvarnames,
 bool apply_mass)  {

  using idx_t = quadgrid_t<vector_t<real_t>>::idx_t;
  real_t N = 0.0, xx = 0.0, yy = 0.0;
  idx_t idx = 0;


  for (std::size_t ivar = 0; ivar < std::size(gvarnames); ++ivar) {
    auto & gvar = vars[getkey(gvarnames, ivar)];
    auto const & dprop = dprops.at (getkey(pvarnames, ivar));

    /*for (idx_t ip = 0; ip <= this->num_particles; ++ip) {
      xx = x[ip];
      yy = y[ip];
      auto icell = grid[ptcl_to_grd[ip]];
      for (idx_t inode = 0; inode < 4; ++inode) {
        N = icell.shp(xx, yy, inode) * dprop[ip];
	gvar[icell.gt(inode)] += N;
      }
    }*/

    p2g_helper_t helper(x.cbegin(), y.cbegin(), M.begin(), gvar.begin(),  ptcl_to_grd.cbegin(), grid.num_rows(), grid.hx(), grid.hy(), dprop.cbegin(), apply_mass);
    #ifdef USE_THRUST
    thrust::counting_iterator<idx_t> first_p(0), last_p(this -> num_particles);
    #else
    range<idx_t> rng (0, this->num_particles);
    range<idx_t>::iterator first_p = rng.begin(), last_p = rng.end();
    #endif

    algorithm_namespace::for_each(first_p, last_p, helper);

  }

  /*if (apply_mass)
    for (std::size_t ivar = 0; ivar < std::size (gvarnames); ++ivar)
      for (idx_t ii = 0; ii < M.size (); ++ii) {
        vars[getkey(gvarnames, ivar)][ii]  /= M[ii];
      }*/
}


template<typename str>
void
particles_t::p2gd
(std::map<std::string, vector_t<real_t>> & vars,
 std::initializer_list<str> const & pxvarnames,
 std::initializer_list<str> const & pyvarnames,
 std::string const &area,
 std::initializer_list<str> const & gvarnames,
 bool apply_mass)  {
  using strlist = std::initializer_list<str> const &;
  p2gd<strlist, strlist>
    (vars, pxvarnames, pyvarnames, area, gvarnames, apply_mass);
}

template<typename GVAR_t, typename PVAR_t, typename P2C_t>
class
p2gd_helper_t{

  using idx_t = particles_t::idx_t;
  const PVAR_t x;
  const PVAR_t y;
  const PVAR_t dpropx;
  const PVAR_t dpropy;
  const PVAR_t dproparea;
  const P2C_t ptcl_to_grd;
  const idx_t nrows;
  const real_t hx;
  const real_t hy;
  const GVAR_t M;
  GVAR_t gvar;

  bool apply_mass;

public :

  p2gd_helper_t (const PVAR_t x_, const PVAR_t y_,  const GVAR_t M_,
		GVAR_t gvar_, const P2C_t ptcl_to_grd_, const idx_t nrows_,
		const real_t hx_, const real_t hy_, const PVAR_t dpropx_, const PVAR_t dpropy_, const PVAR_t dproparea_, bool apply_mass_)
    : x(x_), y(y_), M(M_), gvar(gvar_), 
      ptcl_to_grd(ptcl_to_grd_), nrows(nrows_), hx(hx_), hy(hy_),
      dpropx(dpropx_), dpropy(dpropy_), dproparea(dproparea_), apply_mass(apply_mass_) {};
  
  DEVICE
  void
  operator() (idx_t ip) {
    using qgt = quadgrid_t<GVAR_t>;
    real_t Nx = 0.0, Ny = 0.0;
    auto xx = x[ip];
    auto yy = y[ip];
    auto r = qgt::gind2row (ptcl_to_grd[ip], nrows);
    auto c = qgt::gind2col (ptcl_to_grd[ip], nrows);
    for (idx_t inode = 0; inode < 4; ++inode) {  
      Nx = apply_mass ? qgt::shg (xx, yy, 0, inode, c, r, hx, hy) * M[qgt::gt(inode, c, r, nrows)] :
	    qgt::shg (xx, yy, 0, inode, c, r, hx, hy);
      Ny = apply_mass ? qgt::shg (xx, yy, 1, inode, c, r, hx, hy) * M[qgt::gt(inode, c, r, nrows)] :
      qgt::shg (xx, yy, 1, inode, c, r, hx, hy);

      atomicAdd(&(gvar[qgt::gt(inode, c, r, nrows)]), (Nx*dpropx[ip] + Ny*dpropy[ip])*dproparea[ip]);
    }
  } 


};


template<typename GT, typename PT>
void
particles_t::p2gd
(std::map<std::string, vector_t<real_t>> & vars,
 PT const & pxvarnames,
 PT const & pyvarnames,
 std::string const &area,
 GT const & gvarnames,
 bool apply_mass)  {

  using idx_t = quadgrid_t<vector_t<real_t>>::idx_t;
  real_t xx = 0.0, yy = 0.0, Nx = 0.0, Ny = 0.0;

  for (std::size_t ivar = 0; ivar < std::size (gvarnames); ++ivar) {
    auto & gvar = vars[getkey(gvarnames, ivar)];
    auto const & dpropx = dprops.at (getkey(pxvarnames, ivar));
    auto const & dpropy = dprops.at (getkey(pyvarnames, ivar));
    auto const & dproparea = dprops.at (area);

    p2gd_helper_t helper(x.cbegin(), y.cbegin(), M.cbegin(), gvar.begin(), ptcl_to_grd.cbegin(), grid.num_rows(), grid.hx(), grid.hy(), 
    dpropx.begin(), dpropy.begin(), dproparea.begin(), apply_mass);

    #ifdef USE_THRUST
    thrust::counting_iterator<idx_t> first_p(0), last_p(this -> num_particles);
    #else
    range<idx_t> rng (0, this->num_particles);
    range<idx_t>::iterator first_p = rng.begin(), last_p = rng.end();
    #endif

    algorithm_namespace::for_each(first_p, last_p, helper);

    /*for (idx_t ip = 0; ip <= this->num_particles; ++ip) {
      xx = x[ip];
      yy = y[ip];
      auto icell = grid[ptcl_to_grd[ip]];
      for (idx_t inode = 0; inode < 4; ++inode) {
        Nx = icell.shg (xx, yy, 0, inode);
        Ny = icell.shg (xx, yy, 1, inode);
        gvar[icell.gt(inode)] += (Nx * dpropx[ip] + Ny * dpropy[ip]) * dproparea[ip];
      }
    }*/

  }

  /*if (apply_mass)
    for (std::size_t ivar = 0; ivar < std::size (gvarnames); ++ivar)
      for (idx_t ii = 0; ii < M.size (); ++ii) {
        vars[getkey(gvarnames, ivar)][ii]  /= M[ii];
      }*/

}

template<typename str>
void
particles_t::g2p
(const std::map<std::string, vector_t<real_t>> & vars,
 std::initializer_list<str> const & gvarnames,
 std::initializer_list<str> const & pvarnames,
 bool apply_mass) {
  using strlist = std::initializer_list<str> const &;
  g2p<strlist, strlist> (vars, gvarnames,
                         pvarnames, apply_mass);
}

//! @brief Template class for the implementation
//! of the `g2p` method.
template<typename GVAR_t, typename PVAR_t, typename P2C_t>
class
g2p_helper_t {

  using idx_t = particles_t::idx_t;
  const PVAR_t x;
  const PVAR_t y;
  GVAR_t M;
  const GVAR_t gvar;
  const P2C_t ptcl_to_grd;
  const idx_t nrows;
  const real_t hx;
  const real_t hy;
  PVAR_t dprop;
  bool apply_mass;
  
public :

  g2p_helper_t (const PVAR_t x_, const PVAR_t y_,
		const GVAR_t gvar_, const P2C_t ptcl_to_grd_, const idx_t nrows_,
		const real_t hx_, const real_t hy_, PVAR_t dprop_, bool apply_mass_)
    : x(x_), y(y_), gvar(gvar_),
      ptcl_to_grd(ptcl_to_grd_), nrows(nrows_), hx(hx_), hy(hy_),
      dprop(dprop_), apply_mass(apply_mass_) {};
  
  DEVICE
  void
  operator() (idx_t ip) {
    using qgt = quadgrid_t<GVAR_t>;
    real_t N = 0.0;
    auto xx = x[ip];
    auto yy = y[ip];
    auto r = qgt::gind2row (ptcl_to_grd[ip], nrows);
    auto c = qgt::gind2col (ptcl_to_grd[ip], nrows);
    for (idx_t inode = 0; inode < 4; ++inode) {  
      N = apply_mass ? qgt::shp (xx, yy, inode, c, r, hx, hy) * M[qgt::gt(inode, c, r, nrows)] :
	qgt::shp (xx, yy, inode, c, r, hx, hy);
      dprop[ip] += N * gvar[qgt::gt(inode, c, r, nrows)];
    }
  }  
};


template<typename GT, typename PT>
void
particles_t::g2p
(const std::map<std::string, vector_t<real_t>>& vars,
 GT const & gvarnames,
 PT const & pvarnames,
 bool apply_mass) {

  using idx_t = particles_t::idx_t;
  real_t N = 0.0, xx = 0.0, yy = 0.0;
  idx_t idx = 0;

  for (std::size_t ivar = 0; ivar < std::size (gvarnames); ++ivar) {
    auto & dprop = dprops.at (getkey (pvarnames, ivar));
    auto const & gvar = vars.at (getkey (gvarnames, ivar));
    
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

    g2p_helper_t helper (x.begin (), y.begin (),
    			 gvar.cbegin (), ptcl_to_grd.cbegin (),
    			 grid.num_rows (), grid.hx (), grid.hy (),
    			 dprop.begin (), apply_mass);

    #ifdef USE_THRUST
    thrust::counting_iterator<idx_t> first_p(0), last_p(this -> num_particles);
    #else
    range<idx_t> rng (0, this->num_particles);
    range<idx_t>::iterator first_p = rng.begin(), last_p = rng.end();
    #endif
    
    algorithm_namespace::for_each(first_p, last_p, helper);
    
  }
}

template<typename str>
void
particles_t::g2pd
(const std::map<std::string, vector_t<real_t>>& vars,
 std::initializer_list<str> const & gvarnames,
 std::initializer_list<str> const & pxvarnames,
 std::initializer_list<str> const & pyvarnames,
 bool apply_mass) {
  using strlist = std::initializer_list<str> const &;
  g2pd<strlist, strlist> (vars, gvarnames, pxvarnames,
                          pyvarnames, apply_mass);
}

//! @brief Template class for the implementation
//! `g2pd` method.
template<typename GVAR_t, typename PVAR_t, typename P2C_t>
class
g2pd_helper_t {

  using idx_t = particles_t::idx_t;
  const PVAR_t x;
  const PVAR_t y;
  const GVAR_t M;
  const GVAR_t gvar;
  const P2C_t ptcl_to_grd;
  const idx_t nrows;
  const real_t hx;
  const real_t hy;
  PVAR_t dpropx;
  PVAR_t dpropy;
  bool apply_mass;
  
public :

  g2pd_helper_t (const PVAR_t x_, const PVAR_t y_, const GVAR_t M_,
		 const GVAR_t gvar_, const P2C_t ptcl_to_grd_, const idx_t nrows_,
		 const real_t hx_, const real_t hy_, PVAR_t dpropx_, PVAR_t dpropy_,
		 bool apply_mass_)
    : x(x_), y(y_), gvar(gvar_),
      ptcl_to_grd(ptcl_to_grd_), nrows(nrows_), hx(hx_), hy(hy_),
      dpropx(dpropx_), dpropy(dpropy_), apply_mass(apply_mass_) {};
  
  DEVICE
  void
  operator() (idx_t ip) {
    using qgt = quadgrid_t<GVAR_t>;
    real_t Nx = 0.0, Ny = 0.0;
    auto xx = x[ip];
    auto yy = y[ip];
    auto r = qgt::gind2row (ptcl_to_grd[ip], nrows);
    auto c = qgt::gind2col (ptcl_to_grd[ip], nrows);

    for (idx_t inode = 0; inode < 4; ++inode) {
      Nx = apply_mass ?
	qgt::shg (xx, yy, 0, inode, c, r, hx, hy) * M[qgt::gt(inode, c, r, nrows)] :
	qgt::shg (xx, yy, 0, inode, c, r, hx, hy);
      Ny = apply_mass ?
	qgt::shg (xx, yy, 1, inode, c, r, hx, hy) * M[qgt::gt(inode, c, r, nrows)] :
	qgt::shg (xx, yy, 1, inode, c, r, hx, hy);
      dpropx[ip] += Nx * gvar[qgt::gt(inode, c, r, nrows)];
      dpropy[ip] += Ny * gvar[qgt::gt(inode, c, r, nrows)];

    }
  } 
};
  
template<typename GT, typename PT>
void
particles_t::g2pd
(const std::map<std::string, vector_t<real_t>>& vars,
 GT const & gvarnames,
 PT const & pxvarnames,
 PT const & pyvarnames,
 bool apply_mass) {

  using idx_t = quadgrid_t<vector_t<real_t>>::idx_t;
  real_t Nx = 0.0, Ny = 0.0, xx = 0.0, yy = 0.0;

  for (std::size_t ivar = 0; ivar < std::size (gvarnames); ++ivar) {
    auto const & gvar = vars.at (getkey (gvarnames, ivar));
    auto & dpropx = dprops.at (getkey (pxvarnames, ivar));
    auto & dpropy = dprops.at (getkey (pyvarnames, ivar));

    // for (idx_t ip = 0; ip <= this->num_particles; ++ip) {
    //   xx = x[ip];
    //   yy = y[ip];
    //   auto icell = grid[ptcl_to_grd[ip]];
    //   for (idx_t inode = 0; inode < 4; ++inode) {
    //     Nx = apply_mass ?
    //       icell.shg(xx, yy, 0, inode) * M[icell.gt(inode)] :
    //       icell.shg(xx, yy, 0, inode);
    //     Ny = apply_mass ?
    //       icell.shg(xx, yy, 1, inode) * M[icell.gt(inode)] :
    //       icell.shg(xx, yy, 1, inode);
    //     dpropx[ip] += Nx * gvar[icell.gt(inode)];
    //     dpropy[ip] += Ny * gvar[icell.gt(inode)];
    //   }
    // }

    g2pd_helper_t helper (x.begin (), y.begin (), M.cbegin (),
    			 gvar.cbegin (), ptcl_to_grd.cbegin (),
    			 grid.num_rows (), grid.hx (), grid.hy (),
    			 dpropx.begin (), dpropy.begin (), apply_mass);
     
    #ifdef USE_THRUST
    thrust::counting_iterator<idx_t> first_p(0), last_p(this -> num_particles);
    #else
    range<idx_t> rng (0, this->num_particles);
    range<idx_t>::iterator first_p = rng.begin(), last_p = rng.end();
    #endif
    
    algorithm_namespace::for_each(first_p, last_p, helper);

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
