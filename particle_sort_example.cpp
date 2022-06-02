#include <random>
#include <quadgrid_cpp.h>

using idx_t = quadgrid_t<std::vector<double>>::idx_t;


struct particles_t {

  std::vector<idx_t> label;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> mass;
  std::vector<double> velocity;
  std::vector<double> energy;
  std::vector<double> M;
  std::map<idx_t, std::vector<idx_t>> grd_to_ptcl;
  quadgrid_t<std::vector<double>>& grid;
  
  particles_t (idx_t n, quadgrid_t<std::vector<double>>& grid_)
    : label(n,0), x(n,0.0), y(n,0.0), mass(n,1.0),
      velocity (n,1.0), energy (n,1.0), grid(grid_) {

    M = std::vector<double> (grid.num_global_nodes (), 0.0);
   
     for (auto icell = grid.begin_cell_sweep ();
          icell != grid.end_cell_sweep (); ++icell) {      
       for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode) {
         M[icell->gt (inode)] += (grid.hx () / 2.) * (grid.hy () / 2.);          
       }
     }
  
     std::random_device rd;  
     std::mt19937 gen (rd ()); 
     std::uniform_real_distribution<> dis (0.0, 1.0);
     std::generate (x.begin (), x.end (), [&] () { return dis (gen) * grid.num_cols () * grid.hx (); });
     std::generate (y.begin (), y.end (), [&] () { return dis (gen) * grid.num_rows () * grid.hy (); });
  
     idx_t ilab = 0;
     std::iota (label.begin (), label.end (), ilab);
     
     for (auto ii = 0; ii - x.size (); ++ii) {
       idx_t c = static_cast<idx_t> (std::floor (x[ii] / grid.hx ()));
       idx_t r = static_cast<idx_t> (std::floor (y[ii] / grid.hy ()));
       
       grd_to_ptcl[grid.sub2gind (r, c)].push_back (ii);
     }
  };

  void
  p2g (std::vector<double>& gm, std::vector<double>& gv, std::vector<double>& ge) const {

    for (auto icell = grid.begin_cell_sweep ();
         icell != grid.end_cell_sweep (); ++icell) {
      for (idx_t ii = 0; ii < grd_to_ptcl.at (icell->get_global_cell_idx ()).size (); ++ii) {
        double xx = x[grd_to_ptcl.at(icell->get_global_cell_idx ())[ii]];
        double yy = y[grd_to_ptcl.at(icell->get_global_cell_idx ())[ii]];
        double mm = mass[grd_to_ptcl.at(icell->get_global_cell_idx ())[ii]];
        double vv = velocity[grd_to_ptcl.at(icell->get_global_cell_idx ())[ii]];
        double ee = energy[grd_to_ptcl.at(icell->get_global_cell_idx ())[ii]];
        for (idx_t inode = 0; inode < 4; ++inode) {
          double N = icell->shp(xx, yy, inode);
          gm[icell->t(inode)] += N * mm;
          gv[icell->t(inode)] += N * vv;
          ge[icell->t(inode)] += N * ee;
        }                
      }
    }

    for (idx_t ii = 0; ii < gm.size (); ++ii) {
      gm[ii] /= M[ii];
      gv[ii] /= M[ii];
      ge[ii] /= M[ii];
    }
  };

  void
  g2p (std::vector<double>& mass, std::vector<double>& velocity, std::vector<double>& energy) const {
    
  };

};

int
main (int argc, char *argv[]) {

  quadgrid_t<std::vector<double>> grid;
  grid.set_sizes (4, 4, .25, .25);

  constexpr idx_t num_particles = 1000000;
  particles_t ptcls (num_particles, grid);

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
  
  std::vector<double> rho (grid.num_global_nodes (), 0.0);
  std::vector<double> p (grid.num_global_nodes (), 0.0);
  std::vector<double> ie (grid.num_global_nodes (), 0.0);

  ptcls.p2g (rho, p, ie);
    
  for (auto ii : rho)
    std::cout << ii << std::endl;
  
  return 0;
};


