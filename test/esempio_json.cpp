#include <json.hpp>
#include <particles.h>

#include <fstream>
#include <iostream>
#include <map>
#include <vector>
int
main () {

  using idx_t = quadgrid_t<std::vector<double>>::idx_t;
  quadgrid_t<std::vector<double>> grid;
  grid.set_sizes (4, 6, 1./4., 1./6.);

  constexpr idx_t num_particles = 15;
  particles_t ptcls (num_particles, {"label"}, {"m","vx","vy"}, grid);
  ptcls.dprops["m"].assign (num_particles, (4./16. * 1.) / num_particles );
  ptcls.dprops["vx"].assign (num_particles, 1. );
  ptcls.dprops["vy"].assign (num_particles, 1. );

  idx_t ilabel = 0;
  std::iota (ptcls.iprops["label"].begin (),
	     ptcls.iprops["label"].end (), ilabel);

  {
    std::ofstream jsonfile ("esempio.json");
    jsonfile << std::setw(2) << nlohmann::json (ptcls);
    jsonfile.close ();
  }

  
  nlohmann::json j;
  std::ifstream jsonfile("esempio.json");
  jsonfile >> j;
  quadgrid_t<std::vector<double>> qg (j["grid_properties"]);
  particles_t p (j, qg);
  
  {
    std::ofstream jf ("esempio2.json");
    jf << std::setw(2) << nlohmann::json (p);
    jf.close ();
  }
  
  return 0;
}
