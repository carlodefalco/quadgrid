#include "taylor_dispersion.h"
#include <thrust/extrema.h>

int main(){

#ifndef THRUST_CPU
 int num_gpus;
 auto err = cudaGetDeviceCount (&num_gpus); if (err) return err;
 std::cerr << "num_gpus=" << num_gpus <<std::endl;
 int device;
 //for (int gpu = 0; gpu < num_gpus; ++gpu)
 int gpu = 1;
 {
  int  err = cudaSetDevice (gpu); if (err) return err;
  err = cudaGetDevice (&device); if (err) return err;
  std::cerr << "running on gpu n. " << device << std::endl;}
#endif
 
 using idx_t = particles_t::idx_t;

 cdf::timer::timer_t timer;
 
 // read data from file
 constexpr auto filename = "td.json";
 
 nlohmann::json j;
 std::ifstream inbuf (filename);
 inbuf >> j;

 quadgrid_t<vector_t<real_t>> qg(j["grid_properties"]);
 particles_t p(j,qg);
 
 p.build_mass();
 p.init_particle_mesh();

 std::map<std::string, vector_t<real_t>> vars=
  j["grid_vars"].get<std::map<std::string, vector_t<real_t>>> ();

 p.memcpy_host_to_device();
  
 for (const auto & g : vars)
 p.device_grid_vars[g.first] = g.second;

 inbuf.close ();


 p2g_step1 Jdrift_x(p.device_x.cbegin(), p.device_y.cbegin(), thrust::raw_pointer_cast(p.device_grid_M.data()), 
    thrust::raw_pointer_cast(p.device_grid_vars["Jdrift_x"].data()), p.device_ptcl_to_grd.cbegin(), qg.num_rows(), qg.hx(),
    qg.hy(), p.device_dprops["BETAx"].cbegin(), p.device_dprops["M"].cbegin(), false);

 p2g_step1 Jdrift_y(p.device_x.cbegin(), p.device_y.cbegin(), thrust::raw_pointer_cast(p.device_grid_M.data()),
    thrust::raw_pointer_cast(p.device_grid_vars["Jdrift_y"].data()), p.device_ptcl_to_grd.cbegin(), qg.num_rows(), qg.hx(),
    qg.hy(), p.device_dprops["BETAy"].cbegin(), p.device_dprops["M"].cbegin(), false);

 p2gd_step2 Jdiff_x(p.device_x.cbegin(), p.device_y.cbegin(), thrust::raw_pointer_cast(p.device_grid_M.data()),
   thrust::raw_pointer_cast(p.device_grid_vars["Jdiff_x"].data()), p.device_ptcl_to_grd.cbegin(), qg.num_rows(), qg.hx(),
   qg.hy(), 1.,  p.device_dprops["M"].cbegin(), p.device_dprops["zero"].cbegin(), false);

 p2gd_step2 Jdiff_y(p.device_x.cbegin(), p.device_y.cbegin(), thrust::raw_pointer_cast(p.device_grid_M.data()),
   thrust::raw_pointer_cast(p.device_grid_vars["Jdiff_y"].data()), p.device_ptcl_to_grd.cbegin(), qg.num_rows(), qg.hx(),
   qg.hy(), 1.,  p.device_dprops["zero"].cbegin(), p.device_dprops["M"].cbegin(), false); 

 g2p_step3 VX(thrust::raw_pointer_cast(p.device_x.data()), thrust::raw_pointer_cast(p.device_y.data()), p.device_grid_M.cbegin(), 
    p.device_grid_vars["Jdrift_x"].cbegin(), p.device_grid_vars["Jdiff_x"].cbegin(),  p.device_ptcl_to_grd.cbegin(), qg.num_rows (), qg.hx (), qg.hy (), 
    thrust::raw_pointer_cast(p.device_dprops["VX"].data()), false);

 g2p_step3 VY(thrust::raw_pointer_cast(p.device_x.data()), thrust::raw_pointer_cast(p.device_y.data()), p.device_grid_M.cbegin(), 
    p.device_grid_vars["Jdrift_y"].cbegin(), p.device_grid_vars["Jdiff_y"].cbegin(),  p.device_ptcl_to_grd.cbegin(), qg.num_rows (), qg.hx (), qg.hy (),
    thrust::raw_pointer_cast(p.device_dprops["VY"].data()), false);
 
 boundary bc(thrust::raw_pointer_cast(p.device_grid_vars["Jdiff_y"].data()), qg.num_cols(), qg.num_rows());
 
 stepper step(p.device_x.begin(), p.device_y.begin(), p.device_dprops["VX"].begin(), p.device_dprops["VY"].begin(), 1e-8);

 thrust::counting_iterator<idx_t> first_p(0), last_p(p.num_particles), first_n(0), last_n(qg.num_cols());

 constexpr int nsave = 10;
 constexpr double tmax = 0.02;
 constexpr double dtsave = tmax / nsave;
 double t = 0.;

 // Saving iteration stepping
 for(int isave=0; isave < nsave; ++isave){
  
  //INTERMIDIATE COPY
  p.p2g(p.device_grid_vars, {"M"}, {"rho"}, true);
  p.memcpy_device_to_host();
  //The following copy cannot be included in the memcpy function since vars is not defined in particles.h
  thrust::copy (p.device_grid_vars["rho"].cbegin(), p.device_grid_vars["rho"].cend(), vars.at("rho").begin());
  thrust::copy (p.device_grid_vars["Jdrift_x"].cbegin(), p.device_grid_vars["Jdrift_x"].cend(), vars.at("Jdrift_x").begin());
  thrust::copy (p.device_grid_vars["Jdrift_y"].cbegin(), p.device_grid_vars["Jdrift_y"].cend(), vars.at("Jdrift_y").begin());
  thrust::copy (p.device_grid_vars["Jdiff_x"].cbegin(), p.device_grid_vars["Jdiff_x"].cend(), vars.at("Jdiff_x").begin());
  thrust::copy (p.device_grid_vars["Jdiff_y"].cbegin(), p.device_grid_vars["Jdiff_y"].cend(), vars.at("Jdiff_y").begin());

  // write particle data to file
  const std::string ofilename = "particle";
  const std::string ofileext = ".csv";
  const std::string numfile = std::string(".") + std::to_string(isave);
	
  std::ofstream outbuf (ofilename + numfile + ofileext);  
  p.print<particles_t::output_format::csv> (outbuf);
	
  outbuf.close ();
	
  // write grid data to file
  const std::string gfilename = std::string("grid.") + std::to_string(isave) + std::string(".vts");
  qg.vtk_export (gfilename.c_str(), vars);

  //Time Stepping
 
  while(t < dtsave*(isave +1)){  
  thrust::fill(p.device_grid_vars["Jdrift_x"].begin(), p.device_grid_vars["Jdrift_x"].end(), 0.0);
  thrust::fill(p.device_grid_vars["Jdrift_y"].begin(), p.device_grid_vars["Jdrift_y"].end(), 0.0);
  thrust::fill(p.device_grid_vars["Jdiff_x"].begin(), p.device_grid_vars["Jdiff_x"].end(), 0.0);   
  thrust::fill(p.device_grid_vars["Jdiff_y"].begin(), p.device_grid_vars["Jdiff_y"].end(), 0.0);
  thrust::fill(p.device_grid_vars["rho"].begin(), p.device_grid_vars["rho"].end(), 0.0); 
  thrust::fill(p.device_dprops["BETAx"].begin(), p.device_dprops["BETAx"].end(), 0.0); 
  thrust::fill(p.device_dprops["BETAy"].begin(), p.device_dprops["BETAy"].end(), 0.0); 
  thrust::fill(p.device_dprops["VX"].begin(), p.device_dprops["VX"].end(), 0.0); 
  thrust::fill(p.device_dprops["VY"].begin(), p.device_dprops["VY"].end(), 0.0);  
 
  //G2P of drift velocity
  timer.tic("g2p_1");
  p.g2p(p.device_grid_vars, {"betax", "betay"}, {"BETAx", "BETAy"}, false);
  timer.toc("g2p_1");

  //P2G of Jdrift
  timer.tic("p2g");
  p.p2g(p.device_grid_vars, {"M"}, {"rho"}, false); 
  p.p2g(p.device_grid_vars, Jdrift_x);
  p.p2g(p.device_grid_vars, Jdrift_y);
 
  thrust::transform(p.device_grid_vars["Jdrift_x"].begin(), p.device_grid_vars["Jdrift_x"].end(),
                    p.device_grid_vars["rho"].begin(), p.device_grid_vars["Jdrift_x"].begin(), thrust::divides<double>());

  thrust::transform(p.device_grid_vars["Jdrift_y"].begin(), p.device_grid_vars["Jdrift_y"].end(), 
                    p.device_grid_vars["rho"].begin(), p.device_grid_vars["Jdrift_y"].begin(), thrust::divides<double>());
  timer.toc("p2g");
  
  //P2GD of Jdiff 
  timer.tic("p2gd");
  p.p2gd(p.device_grid_vars, Jdiff_x);
  p.p2gd(p.device_grid_vars, Jdiff_y);

  thrust::transform(p.device_grid_vars["Jdiff_x"].begin(), p.device_grid_vars["Jdiff_x"].end(),
                    p.device_grid_vars["rho"].begin(), p.device_grid_vars["Jdiff_x"].begin(), thrust::divides<double>());

  thrust::transform(p.device_grid_vars["Jdiff_y"].begin(), p.device_grid_vars["Jdiff_y"].end(),
                    p.device_grid_vars["rho"].begin(), p.device_grid_vars["Jdiff_y"].begin(), thrust::divides<double>());
  //set BC
  thrust::for_each(thrust::device, first_n, last_n, bc);
  timer.toc("p2gd");
 
  // G2P of VX and VY
  timer.tic("g2p_2");
  p.g2p(p.device_grid_vars, VX);
  p.g2p(p.device_grid_vars, VY);
  timer.toc("g2p_2");

  //Compute dt
  timer.tic("dt");
  auto maxvx = thrust::max_element (p.device_dprops["VX"].begin (), p.device_dprops["VX"].end ());
  step.dt = .5 * qg.hx() / (*maxvx);
  auto maxvy = thrust::max_element (p.device_dprops["VY"].begin (), p.device_dprops["VY"].end ());
  step.dt = fmin (step.dt, .5 * qg.hy() / (*maxvy));
  if (t + step.dt > dtsave * (isave + 1)) step.dt = dtsave * (isave + 1) - t;
  timer.toc("dt");

  //Moving Particles
  timer.tic("move partcles");
  thrust::for_each(thrust::device, first_p, last_p, step);
  p.update_ptcl_to_grd<particles_t::update_ptcl_to_grd_device>();
  timer.toc("move partcles");
  
  t+=step.dt;

  }
  }
  
  //FINAL COPY   
  p.p2g(p.device_grid_vars, {"M"}, {"rho"}, true);
  p.memcpy_device_to_host();
  //The following copies cannot be included in the memcpy function since vars is not defined in particles.h
  thrust::copy (p.device_grid_vars["rho"].cbegin(), p.device_grid_vars["rho"].cend(), vars.at("rho").begin());
  thrust::copy (p.device_grid_vars["Jdrift_x"].cbegin(), p.device_grid_vars["Jdrift_x"].cend(), vars.at("Jdrift_x").begin());
  thrust::copy (p.device_grid_vars["Jdrift_y"].cbegin(), p.device_grid_vars["Jdrift_y"].cend(), vars.at("Jdrift_y").begin());
  thrust::copy (p.device_grid_vars["Jdiff_x"].cbegin(), p.device_grid_vars["Jdiff_x"].cend(), vars.at("Jdiff_x").begin());
  thrust::copy (p.device_grid_vars["Jdiff_y"].cbegin(), p.device_grid_vars["Jdiff_y"].cend(), vars.at("Jdiff_y").begin());

  // write particle data to file
  const std::string ofilename = "particle";
  const std::string ofileext = ".csv";
  const std::string numfile = std::string(".") + std::to_string(nsave);
	
  std::ofstream outbuf (ofilename + numfile + ofileext);  
  p.print<particles_t::output_format::csv> (outbuf);
	
  outbuf.close ();
	
  // write grid data to file
  const std::string gfilename = std::string("grid.") + std::to_string(nsave) + std::string(".vts");
  qg.vtk_export (gfilename.c_str(), vars);
      

      
 
 timer.print_report();
 return 0; 
}
