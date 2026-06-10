#ifndef CHANNEL_H
#define CHANNEL_H

#include <json.hpp>
#include <particles.h>
 
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <vector>
 
#include "counter.h"
#include <timer.h>
#include <quadgrid_config.h>


template<typename PVAR_t>
class
stepper {
private :
  PVAR_t x;
  PVAR_t y;
  PVAR_t vx;
  PVAR_t vy;

  
public :
  real_t dt; 
  stepper (PVAR_t x_, PVAR_t y_,
	   PVAR_t vx_, PVAR_t vy_,
	   real_t dt_)
    : x(x_), y(y_), vx(vx_), vy(vy_), dt{dt_} { }

  //! @brief call operator applying motion to the n-th particle.
  
  /// overload of the call operator
  /// to apply motion to the n-th particle
  /// use velocity field

  DEVICE
  void operator() (int n) {
 
    //update particles positions
    x[n] += vx[n] * dt; 
    y[n] += vy[n] * dt; 
      
    // Apply boundary conditions (unelastic walls)
      y[n] = fmin (0.1999, fmax (0.001, y[n]));
  } 

};


template<typename PVAR_t, typename GVAR_t, typename P2C_t>
class
p2g_step1{

  using idx_t = particles_t::idx_t;
  const PVAR_t x;
  const PVAR_t y;
  const PVAR_t dprop1;
  const PVAR_t dprop2;
  const P2C_t ptcl_to_grd;
  const idx_t nrows;
  const real_t hx;
  const real_t hy;
  GVAR_t M;
  GVAR_t gvar;
  bool apply_mass;

public:

  p2g_step1 (const PVAR_t x_, const PVAR_t y_, GVAR_t M_, 
		GVAR_t gvar_, const P2C_t ptcl_to_grd_, const idx_t nrows_,
		const real_t hx_, const real_t hy_, const PVAR_t dprop1_, const PVAR_t dprop2_,  bool apply_mass_)
    : x(x_), y(y_), M(M_), gvar(gvar_), 
      ptcl_to_grd(ptcl_to_grd_), nrows(nrows_), hx(hx_), hy(hy_),
      dprop1(dprop1_), dprop2(dprop2_), apply_mass(apply_mass_) {};


  DEVICE
  void
  operator()(idx_t ip){
  
  using qgt = quadgrid_t<GVAR_t>;
    real_t N = 0.0;
    auto xx = x[ip];
    auto yy = y[ip];
    auto r = qgt::gind2row (ptcl_to_grd[ip], nrows);
    auto c = qgt::gind2col (ptcl_to_grd[ip], nrows);
    for (idx_t inode = 0; inode < 4; ++inode) {  
      N = apply_mass ? qgt::shp (xx, yy, inode, c, r, hx, hy)/M[qgt::gt(inode, c, r, nrows)] :
	    qgt::shp (xx, yy, inode, c, r, hx, hy);
      atomicAdd(&(gvar[qgt::gt(inode, c, r, nrows)]), N*dprop1[ip]*dprop2[ip]);
   } 
 }
 
};


template<typename PVAR_t, typename GVAR_t, typename P2C_t>
class
p2gd_step2{

  using idx_t = particles_t::idx_t;
  const PVAR_t x;
  const PVAR_t y;
  const PVAR_t dpropx;
  const PVAR_t dpropy;
  const P2C_t ptcl_to_grd;
  const idx_t nrows;
  const real_t hx;
  const real_t hy;
  const real_t D;
  GVAR_t M;
  GVAR_t gvar;

  bool apply_mass;

public :

  p2gd_step2 (const PVAR_t x_, const PVAR_t y_,  const GVAR_t M_,
		GVAR_t gvar_, const P2C_t ptcl_to_grd_, const idx_t nrows_,
		const real_t hx_, const real_t hy_, const real_t D_, const PVAR_t dpropx_, const PVAR_t dpropy_, bool apply_mass_)
    : x(x_), y(y_), M(M_), gvar(gvar_), 
      ptcl_to_grd(ptcl_to_grd_), nrows(nrows_), hx(hx_), hy(hy_),
      D(D_), dpropx(dpropx_), dpropy(dpropy_),  apply_mass(apply_mass_) {};
  
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
      Nx = apply_mass ? qgt::shg (xx, yy, 0, inode, c, r, hx, hy) / M[qgt::gt(inode, c, r, nrows)] :
	    qgt::shg (xx, yy, 0, inode, c, r, hx, hy);
      Ny = apply_mass ? qgt::shg (xx, yy, 1, inode, c, r, hx, hy) / M[qgt::gt(inode, c, r, nrows)] :
      qgt::shg (xx, yy, 1, inode, c, r, hx, hy);

      atomicAdd(&(gvar[qgt::gt(inode, c, r, nrows)]), (Nx*dpropx[ip] + Ny*dpropy[ip])*D);
    }
  } 

};


template<typename GVAR_t, typename PVAR_t, typename P2C_t>
class
g2p_step3 {

  using idx_t = particles_t::idx_t;
  PVAR_t x;
  PVAR_t y;
  const GVAR_t M;
  const GVAR_t gvar1;
  const GVAR_t gvar2;
  const P2C_t ptcl_to_grd;
  const idx_t nrows;
  const real_t hx;
  const real_t hy;
  PVAR_t dprop;
  bool apply_mass;
  
public :

  g2p_step3 (PVAR_t x_,  PVAR_t y_, const GVAR_t M_,
		const GVAR_t gvar1_, const GVAR_t gvar2_,  const P2C_t ptcl_to_grd_, const idx_t nrows_,
		const real_t hx_, const real_t hy_, PVAR_t dprop_, bool apply_mass_)
    : x(x_), y(y_), M(M_), gvar1(gvar1_), gvar2(gvar2_),
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
      dprop[ip] += N * (gvar1[qgt::gt(inode, c, r, nrows)] + gvar2[qgt::gt(inode, c, r, nrows)]); 
    }
  }  
};


template<typename GVAR_t, typename PVAR_t, typename P2C_t>
class
g2pd_step4 {

  using idx_t = particles_t::idx_t;
  PVAR_t x;
  PVAR_t y;
  const GVAR_t M;
  const GVAR_t gvar1, gvar2, gvar3, gvar4;
  const P2C_t ptcl_to_grd;
  const idx_t nrows;
  const real_t hx;
  const real_t hy;
  PVAR_t dprop;
  bool apply_mass;
  
public :

  g2pd_step4 (PVAR_t x_,  PVAR_t y_, const GVAR_t M_,
		 const GVAR_t gvar1_, const GVAR_t gvar2_, const GVAR_t gvar3_, const GVAR_t gvar4_,  const P2C_t ptcl_to_grd_, const idx_t nrows_,
		 const real_t hx_, const real_t hy_, PVAR_t dprop_, 
		 bool apply_mass_)
    : x(x_), y(y_), M(M_), gvar1(gvar1_), gvar2(gvar2_), gvar3(gvar3_), gvar4(gvar4_),
      ptcl_to_grd(ptcl_to_grd_), nrows(nrows_), hx(hx_), hy(hy_),
      dprop(dprop_), apply_mass(apply_mass_) {};
  
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
      dprop[ip] += Nx * (gvar1[qgt::gt(inode, c, r, nrows)] + gvar2[qgt::gt(inode, c, r, nrows)]) + 
                    Ny * (gvar3[qgt::gt(inode, c, r, nrows)] + gvar4[qgt::gt(inode, c, r, nrows)]);;

    }
  } 
};

template<typename PVAR_t>
class
updateRho{
 PVAR_t rho, divV;
 real_t dt; 
 using idx_t = particles_t::idx_t;
public:
 updateRho(PVAR_t rho_, PVAR_t divV_, real_t dt_):
       rho(rho_), divV(divV_), dt(dt_){};
 
 DEVICE
 void
 operator() (idx_t ip){
    rho[ip] = rho[ip]/(1 + dt*divV[ip]);
 }

};

template<typename GVAR_t>
class boundary{
     using idx_t = particles_t::idx_t;
     GVAR_t Jdiffy;
     idx_t nx, ny;
public:

   boundary(GVAR_t Jdiffy_, idx_t nx_, idx_t ny_) :
     Jdiffy(Jdiffy_), nx(nx_), ny(ny_)  {}
   
   DEVICE
   void
   operator()(idx_t i){
      using qgt = quadgrid_t<GVAR_t>;

      Jdiffy[qgt::sub2gind(ny, i, ny+1)] = 0; 
      Jdiffy[qgt::sub2gind(0, i, ny+1)] = 0;
  
 }   

};



#endif

