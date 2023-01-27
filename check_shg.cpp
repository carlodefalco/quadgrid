#include <algorithm>
#include <fstream>
#include <particles.h>
#include <quadgrid_cpp.h>


using idx_t = quadgrid_t<std::vector<double>>::idx_t;

int
main (int argc, char *argv[]) {

  quadgrid_t<std::vector<double>> grid;
  constexpr double hx = 1.;
  constexpr double hy = 2.;
  
  grid.set_sizes (1, 1, hx, hy);

  constexpr idx_t num_particles = 100;

  std::vector<double>
    x(num_particles);
  std::vector<double>
    y(num_particles);
  
  {
    idx_t k = 0;
    for(idx_t i=0; i<10; ++i) {
      for(idx_t j=0; j<10; ++j) {
	x[k] = i*hx/10 + hx/10/2;
	y[k++] = j*hy/10 + hy/10/2;
      }
    }
  }

  
  particles_t ptcls (num_particles, {"label"}, {"shp", "shgx", "shgy"}, grid, x, y);
  
  
  ptcls.dprops["shp"].resize (num_particles);
  ptcls.dprops["shgx"].resize (num_particles);
  ptcls.dprops["shgy"].resize (num_particles);
  
  ptcls.init_particle_mesh ();
  ptcls.build_mass ();

  for (auto icell = grid.begin_cell_sweep ();
       icell != grid.end_cell_sweep (); ++icell) {
    
    for (idx_t inode = 0; inode < 4; ++inode) {
      for (idx_t i = 0; i < num_particles; ++i) {
	ptcls.dprops["shp"][i] = icell->shp(ptcls.x[i], ptcls.y[i], inode);
	ptcls.dprops["shgx"][i] = icell->shg(ptcls.x[i], ptcls.y[i], 0, inode);
	ptcls.dprops["shgy"][i] = icell->shg(ptcls.x[i], ptcls.y[i], 1, inode);
      }
      std::string filename = "shp.csv";
      filename = std::to_string(inode) + filename;
      std::ofstream of(filename.c_str ());
      ptcls.print<particles_t::output_format::csv> (of);
      of.close ();
    }
    
  }


  return 0;
};
