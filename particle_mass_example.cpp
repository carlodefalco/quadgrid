#include <algorithm>
#include <random>
#include <quadgrid_cpp.h>
#include <particles.h>


using idx_t = quadgrid_t<std::vector<double>>::idx_t;

int
main (int argc, char *argv[]) {

  quadgrid_t<std::vector<double>> grid;
  grid.set_sizes (40, 160, 1./160., 1./160.);

  constexpr idx_t num_particles = 1000000;
  particles_t ptcls (num_particles, {"label"}, {"m","vx","vy"}, grid);
  std::cerr << " 1 " << "\n";

  ptcls.dprops["m"].assign (num_particles, (4./16. * 1.) / num_particles );
  std::cerr << " 2 " << "\n";


  idx_t ilabel = 0;
  std::iota (ptcls.iprops["label"].begin (), ptcls.iprops["label"].end (), ilabel);
  std::cerr << " 3 " << "\n";


  // for (auto ii : ptcls.x)
  //    std::cout << ii << " ";
  // std::cout << std::endl << std::endl;
  
  
  //  for (auto icell = grid.begin_cell_sweep ();
  //       icell != grid.end_cell_sweep (); ++icell) {
  //    if (ptcls.grd_to_ptcl.count (icell->get_global_cell_idx ())) {             
  //      auto & plist = ptcls.grd_to_ptcl.at (icell->get_global_cell_idx ());
  //      for (auto ii = 0; ii < plist.size (); ++ii) {
  // 	auto jj = plist[ii];
  // 	if (ptcls.x[jj] > .5) {
  // 	  ptcls.x[jj] = ptcls.x[jj] - .5;
  // 	}
  //    }
  //   }
  //  }
  //  for (auto ii : ptcls.x)
  //    std::cout << ii << " ";
  //  std::cout << std::endl << std::endl;


  
  // for (auto & icell : ptcls.grd_to_ptcl) {
  //   auto cellnum = icell.first;
  //   auto & plist = icell.second;
  //   for (auto jj : plist) {
  //     if (ptcls.x[jj] > .5) {
  // 	//std::cerr << " jj = " << jj << " ptcls.x[jj] = " << ptcls.x[jj];
  // 	ptcls.dprops["m"][jj] = .0;
  // 	//std::cerr << " ptcls.x[jj] = " << ptcls.x[jj] << std::endl;
  //     }      
  //   }
  // }

  for (idx_t ip = 0; ip < ptcls.x.size (); ++ip) {    
    if (ptcls.x[ip] > .5) {
      ptcls.dprops["m"][ip] = .0;
    }      
  }


  //  ptcls.init_particle_mesh ();

  std::cerr << " 4 " << "\n";

  std::map<std::string, std::vector<double>>
    vars{{"m", std::vector<double>(grid.num_global_nodes (), 0.0)}};
  std::cerr << " 5 " << "\n";

  std::vector<std::string> pvarnames{"m"};
  ptcls.p2g (vars, pvarnames, pvarnames, true);
  std::cerr << " 6 " << "\n";
  
  //  for (auto ii : vars["m"])
  //  std::cout << ii << std::endl;

  grid.vtk_export ("particles_mass_example.vts", vars);
  
  std::cerr << " 7 " << "\n";
  
  return 0;
};


