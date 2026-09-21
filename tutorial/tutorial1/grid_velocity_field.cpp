#include <json.hpp>
#include <particles.h>

#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#include "counter.h"
#include <timer.h>


/// \page tutorial Tutorial
/// This is a tutorial, see the annotated source 
/// source at .grid_velocity_field.cpp

//! @brief Functor class for moving particles.

/// This class captures references to particle
/// positions and velocities and overloads the
/// call operator to allow use in STL algorithms
/// in particular the std::for_each method 
class
stepper {
private :
  std::vector<double> &x;
  std::vector<double> &y;
  std::vector<double> &vx;
  std::vector<double> &vy;

  double dt, D; 
  std::function<double (void)> normal;
  
public :
  
  stepper (std::vector<double> &x_, std::vector<double> &y_,
	   std::vector<double> &vx_, std::vector<double> &vy_,
	   double dt_, double D_, std::function<double (void)> &noise_)
    : x(x_), y(y_), vx(vx_), vy(vy_), dt{dt_}, D{D_}, normal(noise_) { }

  //! @brief call operator applying motion to the n-th particle.
  
  /// overload of the call operator
  /// to apply motion to the n-th particle
  /// use velocity field for deterministic
  /// motion component plus gaussian noise
  /// to represent diffusion/Brownian motion
  void operator() (int n) {
    double dxb, dyb;
  
    //Brownian motion displacements
    dxb=std::sqrt (2*D*dt) * normal();
    dyb=std::sqrt (2*D*dt) * normal();
 
    //update particles positions
    x[n] += vx[n] * dt + dxb;
    y[n] += vy[n] * dt + dyb;
      
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
  inbuf >> j;

  // create grid from fields in a json object
  quadgrid_t<std::vector<double>> qg (j["grid_properties"]);

  // create particles from properties in the json object
  // and the above created grid
  particles_t p (j, qg);
  
  // as we will use the mass matrix we must initialize it manually
  p.build_mass ();

  // the variables defined on the grid are not class members
  std::map<std::string, std::vector<double>> vars= 
    j["grid_vars"].get<std::map<std::string, std::vector<double>>> ();
    
  inbuf.close ();

  // Diffusion is modelled as a Gaussian process
  std::random_device rd2;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen2;
  std::normal_distribution<> normal;
  std::function<double ()> noise = [&gen2, &rd2, &normal] () { return normal(gen2); };

  // Create the callable object to be used for moving the particles
  // capture references to the particle positions and velocities
  stepper state (p.x, p.y, p.dprops["VX"], p.dprops["VY"], 1., 1.e-5, noise);

  // Create particle <-> grid connectivity
  // must be updated explicitely
  p.init_particle_mesh ();


  // Time stepping
  for (int it = 0; it < 1000; ++it) {

    // Clean up grid variables at each step!
    std::fill(p.dprops["VX"].begin (), p.dprops["VX"].end (), 0.0);
    std::fill(p.dprops["VY"].begin (), p.dprops["VY"].end (), 0.0);
    std::fill(vars["rho"].begin (), vars["rho"].end (), 0.0);

    // G2P : interpolate velocity at particle positions
    timer.tic("g2p");
    p.g2p (vars, {"vx", "vy"}, {"VX", "VY"});
    timer.toc("g2p");

    // Move particles
    timer.tic("move partcles");

    // You can use a loop
    // for (int ip = 0; ip < p.num_particles; ++ip) {
    //  state(ip);
    // }

    // Or use an STL algorithm
    range rng (0, p.num_particles);
    std::for_each (rng.begin (), rng.end (), state);
      
    timer.toc("move partcles");

    // Rebuild particle <-> grid connectivity
    // must be updated explicitely
    timer.tic("init_particle_mesh");
    p.init_particle_mesh ();
    timer.toc("init_particle_mesh");
    
    // Project particle masses onto the greed and
    // build a density field, only used for output
    timer.tic("p2g");
    p.p2g (vars, {"M"}, {"rho"}, true);
    timer.toc("p2g");

    // Don't save at every timestep
    if (it % 50 == 0) {
      
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
      
    }
  }

  // print timing information
  timer.print_report();
  return 0;
}
