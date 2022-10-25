#include <algorithm>
#include <random>
#include <quadgrid_cpp.h>
#include <particles.h>
#include <map>
#include <iostream>

using idx_t = quadgrid_t<std::vector<double>>::idx_t;

int
main (int argc, char *argv[]) {

  quadgrid_t<std::vector<double>> grid;
  grid.set_sizes (16, 16, 1./16., 1./16.);

  constexpr idx_t num_particles = 1000000;
  particles_t ptcls (num_particles, {"label"}, {"m", "vx", "vy"}, grid);
  ptcls.dprops["m"].assign (num_particles, 1. / static_cast<double>(num_particles));

  idx_t ilabel = 0;
  std::iota (ptcls.iprops["label"].begin (), ptcls.iprops["label"].end (), ilabel);
  /*
  //
  // This will produce very verbose output
  //
  for (auto icell = grid.begin_cell_sweep ();
       icell != grid.end_cell_sweep (); ++icell) {

    std::cout << "cell " << icell->get_global_cell_idx ()
              << " xlims [" << icell->p (0, 0) << ", "
              << icell->p (0, 3) << "]" 
              << " ylims [" << icell->p (1, 0) << ", "
              << icell->p (1, 3) << "]" << std::endl;

    for (auto ii = 0; ii < ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++ii) {
      std::cout << "\tparticle " << ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[ii]
                << ": (" << ptcls.x[ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[ii]]
                << ", " << ptcls.y[ptcls.grd_to_ptcl.at(icell->get_global_cell_idx ())[ii]]
                << ")" << std::endl;
    }

  }
  */

  std::map<std::string, std::vector<double>>
    vars{{"m", std::vector<double>(grid.num_global_nodes (), 0.0)},
      {"vx", std::vector<double>(grid.num_global_nodes (), 0.0)},
        {"vy", std::vector<double>(grid.num_global_nodes (), 0.0)}};

  ptcls.p2g (vars);
    
  for (auto ii : vars["m"])
    std::cout << ii << std::endl;
  
  return 0;
};


