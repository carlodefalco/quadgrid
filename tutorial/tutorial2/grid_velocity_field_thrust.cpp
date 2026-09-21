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

/// \page tutorial Tutorial
/// This is a tutorial, see the annotated source 
/// source at .grid_velocity_field.cpp

//! @brief Functor class for moving particles.

/// This class captures references to particle
/// positions and velocities and overloads the
/// call operator to allow use in STL algorithms
/// in particular the std::for_each method 
template<typename PVAR_t>
class
stepper {
private :
  PVAR_t x;
  PVAR_t y;
  PVAR_t vx;
  PVAR_t vy;

  real_t dt, D; 
  //std::function<real_t (void)> normal;  DOESNT WORK ON DEVICE
  
public :
  
  stepper (PVAR_t x_, PVAR_t y_,
	   PVAR_t vx_, PVAR_t vy_,
	   real_t dt_, real_t D_)
    : x(x_), y(y_), vx(vx_), vy(vy_), dt{dt_}, D{D_} { }

  //! @brief call operator applying motion to the n-th particle.
  
  /// overload of the call operator
  /// to apply motion to the n-th particle
  /// use velocity field for deterministic
  /// motion component plus gaussian noise
  /// to represent diffusion/Brownian motion
  DEVICE
  void operator() (int n) {
    //real_t dxb, dyb;
  
    //Brownian motion displacements
    //dxb=std::sqrt (2*D*dt) * normal();
    //dyb=std::sqrt (2*D*dt) * normal();
 
    //update particles positions
    x[n] += vx[n] * dt; //+ dxb;
    y[n] += vy[n] * dt; //+ dyb;
      
    // Apply boundary conditions (unelastic walls)
    y[n] = std::min (1.999, std::max (0.001, y[n]));
    x[n] = std::min (1.999, std::max (0.001, x[n]));
  } 

};


//! @brief main implementing the time loop.

int
main () {

  cdf::timer::timer_t timer;
    
  // read data from file
  constexpr auto filename = "velocity.json";

  nlohmann::json j;
  std::ifstream inbuf (filename);

// Declaration of host and device variables (only if thrust is used)
  #ifdef USE_THRUST
  inbuf >> j;
  using idx_t = particles_t::idx_t;
  std::unique_ptr<quadgrid_t<vector_t<real_t>>> qg;
  std::unique_ptr<particles_t> p;
  qg = std::make_unique<quadgrid_t<vector_t<real_t>>> (j["grid_properties"]);
  p = std::make_unique<particles_t> (j, *qg);


  std::map<std::string, vector_t<real_t>> vars=
    j["grid_vars"].get<std::map<std::string, vector_t<real_t>>> ();

  //Copy from host to device
  p->memcpy_host_to_device();
  p->build_mass();
  p->init_particle_mesh();

  for (const auto & g : vars)
    p->device_grid_vars[g.first] = g.second;

  inbuf.close ();
  
  // Create the callable object to be used for moving the particles
  // capture references to the particle positions and velocities
  stepper state (thrust::raw_pointer_cast(p->device_x.data()), thrust::raw_pointer_cast(p->device_y.data()), thrust::raw_pointer_cast(p->device_dprops["VX"].data()), thrust::raw_pointer_cast(p->device_dprops["VY"].data()), 1., 1.e-5);

// Create iterator for thrust::for_each   
  thrust::counting_iterator<idx_t> first_p(0), last_p(p -> num_particles);

  #else  
  inbuf >> j;

  // create grid from fields in a json object
  quadgrid_t<vector_t<real_t>> qg (j["grid_properties"]);

  // create particles from properties in the json object
  // and the above created grid
  particles_t p (j, qg);
  
  // as we will use the mass matrix we must initialize it manually
  p.build_mass ();

  // the variables defined on the grid are not class members
  std::map<std::string, vector_t<real_t>> vars= 
    j["grid_vars"].get<std::map<std::string, vector_t<real_t>>> ();
    
  inbuf.close ();

  // Diffusion is modelled as a Gaussian process
  /*std::random_device rd2;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen2;
  std::normal_distribution<> normal;
  std::function<real_t ()> noise = [&gen2, &rd2, &normal] () { return normal(gen2); };*/

  // Create the callable object to be used for moving the particles
  // capture references to the particle positions and velocities
  stepper state (p.x.data(), p.y.data(), p.dprops["VX"].data(), p.dprops["VY"].data(), 1., 1.e-5);

  // Create particle <-> grid connectivity
  // must be updated explicitely
  p.init_particle_mesh ();

  #endif


  // Time stepping
  for (int it = 0; it < 1000; ++it) {

    #ifdef USE_THRUST
    thrust::fill(p->device_dprops["VX"].begin (), p->device_dprops["VX"].end (), 0.0);
    thrust::fill(p->device_dprops["VX"].begin (), p->device_dprops["VY"].end (), 0.0);
    thrust::fill(p->device_grid_vars["rho"].begin (), p->device_grid_vars["rho"].end (), 0.0);
    #else
    // Clean up grid variables at each step!
    std::fill(p.dprops["VX"].begin (), p.dprops["VX"].end (), 0.0);
    std::fill(p.dprops["VY"].begin (), p.dprops["VY"].end (), 0.0);
    std::fill(vars["rho"].begin (), vars["rho"].end (), 0.0);
    #endif

    // G2P : interpolate velocity at particle positions
    timer.tic("g2p");
    #ifdef USE_THRUST
    p->g2p(p->device_grid_vars, {"vx", "vy"}, {"VX", "VY"});
    #else
    p.g2p (vars, {"vx", "vy"}, {"VX", "VY"});
    #endif

    timer.toc("g2p");

    // Move particles
    timer.tic("move partcles");

    // You can use a loop
    // for (int ip = 0; ip < p.num_particles; ++ip) {
    //  state(ip);
    // }

    // Or use an STL algorithm
    #ifdef USE_THRUST
    thrust::for_each(thrust::device, first_p, last_p, state);
    #else
    // Or use an STL algorithm
    range rng (0, p.num_particles);
    std::for_each (rng.begin (), rng.end (), state);
    #endif

      
    timer.toc("move partcles");

    // Rebuild particle <-> grid connectivity
    // must be updated explicitely
    timer.tic("init_particle_mesh");
    #ifdef USE_THRUST
    p->update_ptcl_to_grd<particles_t::update_ptcl_to_grd_device>();
    #else
    p.update_ptcl_to_grd<particles_t::update_ptcl_to_grd_host>();
    #endif

    timer.toc("init_particle_mesh");
    
    // Project particle masses onto the greed and
    // build a density field, only used for output
    timer.tic("p2g");
    #ifdef USE_THRUST
    p->p2g (p->device_grid_vars, {"M"}, {"rho"}, true);
    #else
    p.p2g (vars, {"M"}, {"rho"}, true);
    #endif

    timer.toc("p2g");

    // Don't save at every timestep
    if (it % 50 == 0) {
      #ifdef USE_THRUST  
      p->memcpy_device_to_host();
      //The following copy cannot be included in the memcpy function since vars is not defined in particles.h
      thrust::copy (p->device_grid_vars["rho"].cbegin(), p->device_grid_vars["rho"].cend(), vars.at("rho").begin());

      // write particle data to file
      const std::string ofilename = "particle";
      const std::string ofileext = ".csv";
      const std::string numfile = std::string(".") + std::to_string(it);
	
      std::ofstream outbuf (ofilename + numfile + ofileext);  
      p->print<particles_t::output_format::csv> (outbuf);
	
      outbuf.close ();
	
      // write grid data to file
      const std::string gfilename = std::string("grid.") + std::to_string(it) + std::string(".vts");
      qg->vtk_export (gfilename.c_str(), vars);
      #else      
      // write particle data to file
      const std::string ofilename = "particle";
      const std::string ofileext = ".csv";
      const std::string numfile = std::string(".") + std::to_string(it);
	
      std::ofstream outbuf (ofilename + numfile + ofileext);  
      p.print<particles_t::output_format::csv> (outbuf);
	
      outbuf.close ();
	
      // write grid data to file
      const std::string gfilename = std::string("grid.") + std::to_string(it) + std::string(".vts");
      qg.vtk_export (gfilename.c_str(), vars);
      #endif
      
    }
  }

  // print timing information
  timer.print_report();
  return 0;
}
