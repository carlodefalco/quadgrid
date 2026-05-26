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
/// source at .grid_velocity_field_comsol.cpp

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
  std::vector<double> &dist;
  std::vector<double> &dist_dx;
  std::vector<double> &dist_dy;
  double hx, hy;
  
  double D; 
  std::function<double (void)> normal;
  
public :

  double dt;
  stepper (std::vector<double> &x_, std::vector<double> &y_,
	   std::vector<double> &vx_, std::vector<double> &vy_,
	   std::vector<double> &dist_,  std::vector<double> &dist_dx_,
	   std::vector<double> &dist_dy_, const double hx_, const double hy_,
	   double dt_, double D_, std::function<double (void)> &noise_)
    : x(x_), y(y_), vx(vx_), vy(vy_), dist(dist_), dist_dx(dist_dx_),
      dist_dy(dist_dy_), hx{hx_}, hy{hy_}, dt{dt_}, D{D_}, normal(noise_) { }

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
    double dx = vx[n] * dt + dxb;
    double dy = vy[n] * dt + dyb;
    // Apply boundary conditions 
    if (dist[n] <= 2*dx || dist[n] < 2*dy) {
      double normdist = std::sqrt (dist_dx[n]*dist_dx[n] + dist_dy[n]*dist_dy[n]);
      double verx = dist_dx[n] / normdist;
      double very = dist_dy[n] / normdist;
      if (dist[n] < 0) { verx = x[n] > 1.6 ? -1. : 1.; very = 0.; }       
      dx = dist_dx[n]/normdist * std::max(0., dist_dx[n]*dx/normdist);
      dy = dist_dy[n]/normdist * std::max(0., dist_dy[n]*dy/normdist);
    }
    
    x[n] += dx;    
    y[n] += dy;
    

    // Constrain particles within domain
    if (x[n] > 8.9) { x[n] = .1; y[n] = 1.6; };
    y[n] = std::min (3.1, std::max (0.1, y[n]));
    x[n] = std::min (8.9, std::max (0.1, x[n]));
  } 

};


//! @brief main implementing the time loop.

int
main () {

  cdf::timer::timer_t timer;
    
  // read data from file
  constexpr auto filename = "velocity_comsol.json";

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


  // Create particle <-> grid connectivity
  // must be updated explicitely
  p.init_particle_mesh ();


  // Create the callable object to be used for moving the particles
  // capture references to the particle positions and velocities
  stepper state (p.x, p.y, p.dprops["VX"], p.dprops["VY"], p.dprops["DIST"],
		 p.dprops["DIST_DX"], p.dprops["DIST_DY"], qg.hx(), qg.hy(),
		 1.e-3, 1.e-1, noise);

  p.g2p  (vars, {"vx", "vy", "dist"}, {"VX", "VY", "DIST"});
  
  constexpr int nsave = 300;
  constexpr double tmax = 2 * 8.5 / 40.;
  constexpr double dtsave = tmax / nsave;
  double t = 0.;
  
  // Time stepping
  for (int isave = 0; isave < nsave; ++isave) {

    // Project particle masses onto the greed and
    // build a density field, only used for output
    timer.tic("p2g");
    p.p2g (vars, {"M"}, {"rho"}, true);
    timer.toc("p2g");

    timer.tic("save particles");
    // write particle data to file
    const std::string ofilename = "particle";
    const std::string ofileext = ".csv";
    const std::string numfile = std::string(".") + std::to_string(isave);
	
    std::ofstream outbuf (ofilename + numfile + ofileext);  
    p.print<particles_t::output_format::csv> (outbuf);	
    outbuf.close ();
    timer.toc("save particles");
    
    // write grid data to file
    timer.tic("save grid");
    const std::string gfilename = std::string("grid.") + std::to_string(isave) + std::string(".vts");
    qg.vtk_export (gfilename.c_str(), vars);
    timer.toc("save grid");

    while (t < dtsave * (isave + 1)) {
      
      // Clean up grid variables at each step!

      timer.tic("clean particle properties interpolated from grid");
      std::fill(p.dprops["VX"].begin (), p.dprops["VX"].end (), 0.0);
      std::fill(p.dprops["VY"].begin (), p.dprops["VY"].end (), 0.0);
      std::fill(p.dprops["DIST"].begin (), p.dprops["DIST"].end (), 0.0);
      std::fill(p.dprops["DIST_DX"].begin (), p.dprops["DIST_DX"].end (), 0.0);
      std::fill(p.dprops["DIST_DY"].begin (), p.dprops["DIST_DY"].end (), 0.0);
      timer.toc("clean particle properties interpolated from grid");

      timer.tic("clean rho");
      std::fill(vars["rho"].begin (), vars["rho"].end (), 0.0);
      timer.toc("clean rho");
    
      // G2P : interpolate velocity and sgd at particle positions
      timer.tic("g2p");
      p.g2p  (vars, {"vx", "vy", "dist"}, {"VX", "VY", "DIST"});
      timer.toc("g2p");

      // Compute dt
      timer.tic("dt");
      auto maxvx = std::max_element (p.dprops["VX"].begin (), p.dprops["VX"].end ());
      state.dt = .5 * qg.hx() / (*maxvx);
      auto maxvy = std::max_element (p.dprops["VY"].begin (), p.dprops["VY"].end ());
      state.dt = std::min (state.dt, .5 * qg.hy() / (*maxvy));
      if (t + state.dt > dtsave * (isave + 1)) state.dt = dtsave * (isave + 1) - t;
      std::cout << "hx = " << qg.hx() << " hy = " << qg.hy() << " vx = " << (*maxvx) << " vy = " << (*maxvy) << std::endl;
      std::cout << "t = " << t << " dt = " << state.dt << " t + dt = " << t + state.dt << std::endl;
      timer.toc("dt");
      
      // G2PD : compute direction away from boundary
      timer.tic("g2pd");
      p.g2pd (vars, {"dist"}, {"DIST_DX"}, {"DIST_DY"});
      timer.toc("g2pd");

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

      t += state.dt; 
    }
  }

  // Project particle masses onto the greed and
  // build a density field, only used for output
  timer.tic("p2g");
  p.p2g (vars, {"M"}, {"rho"}, true);
  timer.toc("p2g");
    
  timer.tic("save particles");
  // write particle data to file
  const std::string ofilename = "particle";
  const std::string ofileext = ".csv";
  const std::string numfile = std::string(".") + std::to_string(nsave);
	
  std::ofstream outbuf (ofilename + numfile + ofileext);  
  p.print<particles_t::output_format::csv> (outbuf);	
  outbuf.close ();
  timer.toc("save particles");
    
  // write grid data to file
  timer.tic("save grid");
  const std::string gfilename = std::string("grid.") + std::to_string(nsave) + std::string(".vts");
  qg.vtk_export (gfilename.c_str(), vars);
  timer.toc("save grid");
    
  // print timing information
  timer.print_report();
  return 0;
}
