#ifndef PARTICLES_H
#define PARTICLES_H

#include <quadgrid_cpp.h>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <random>
#include <unordered_map>

struct
particles_t {

  using idx_t = quadgrid_t<std::vector<double>>::idx_t;
  std::vector<double> x;
  std::vector<double> y;

  std::map<std::string, std::vector<idx_t>> iprops;
  std::map<std::string, std::vector<double>> dprops;

  std::vector<double> M;
  std::map<idx_t, std::vector<idx_t>> grd_to_ptcl;
  const quadgrid_t<std::vector<double>>& grid;

  particles_t (idx_t n, const std::vector<std::string>& ipropnames,
               const std::vector<std::string>& dpropnames,
               const quadgrid_t<std::vector<double>>& grid_)
    : x(n, 0.0), y(n, 0.0), grid(grid_) {

    for (idx_t ii = 0; ii < ipropnames.size (); ++ii) {
      iprops[ipropnames[ii]].assign (n, 0);
    }

    for (idx_t ii = 0; ii < dpropnames.size (); ++ii) {
      dprops[dpropnames[ii]].assign (n, 0.0);
    }

    M = std::vector<double> (grid.num_global_nodes (), 0.0);
    build_mass ();

    random_particle_positions (n);

    init_particle_mesh ();
  };

  void
  init_particle_mesh () {
    for (auto ii = 0; ii - x.size (); ++ii) {
      idx_t c = static_cast<idx_t> (std::floor (x[ii] / grid.hx ()));
      idx_t r = static_cast<idx_t> (std::floor (y[ii] / grid.hy ()));

      grd_to_ptcl[grid.sub2gind (r, c)].push_back (ii);
    }
  };

  void
  random_particle_positions (idx_t n) {
    std::random_device rd;
    std::mt19937 gen (rd ());
    std::uniform_real_distribution<> dis (0.0, 1.0);
    std::generate (x.begin (), x.end (),
                   [&] () { return dis (gen) * grid.num_cols () * grid.hx (); });
    std::generate (y.begin (), y.end (),
                   [&] () { return dis (gen) * grid.num_rows () * grid.hy (); });
  };

  void
  build_mass () {
    M.assign (M.size (), 0.0);
    for (auto icell = grid.begin_cell_sweep ();
         icell != grid.end_cell_sweep (); ++icell) {
      for (auto inode = 0;
           inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell;
           ++inode) {
        M[icell->gt (inode)] += (grid.hx () / 2.) * (grid.hy () / 2.);
      }
    }
  };


	void 
	transfer (const std::string & transfer_type,
						const std::string & variable_name,
						std::map<std::string, std::vector<double>>& vars, //auxiliary grid vars
						const std::vector<double> & input_var, //could be a combination of grid XOR ptcls vars
						bool apply_mass = false) {
		double N = 0.0, xx = 0.0, yy = 0.0;
    idx_t idx = 0;
		const static std::unordered_map<std::string,int> um =
		{
			{"p2g", 1},
			{"p2g_dx", 2},
			{"p2g_dy", 3},
			{"g2p", 4},
			{"g2p_dx", 5},
			{"g2p_dy", 6},
		};
		int transfer_code = um.at(transfer_type);
    for (auto icell = grid.begin_cell_sweep ();
         icell != grid.end_cell_sweep (); ++icell) {
      if (grd_to_ptcl.count (icell->get_global_cell_idx ()) > 0)
        for (idx_t ii = 0;
             ii < grd_to_ptcl.at (icell->get_global_cell_idx ()).size ();
             ++ii) {
          idx = grd_to_ptcl.at(icell->get_global_cell_idx ())[ii];
          xx = x[idx];
          yy = y[idx];
          for (idx_t inode = 0; inode < 4; ++inode) {
            switch (transfer_code) {
							case 1: //p2g
								N = icell->shp(xx, yy, inode);
								vars.at(variable_name).at(icell->t(inode)) += N * input_var.at(ii);
							break;
							case 2: //p2g dx
								N = icell->shg(xx,yy,0,inode);
								vars.at(variable_name).at(icell->t(inode)) += N * input_var.at(ii);
							break;
							case 3: //p2g dy
								N = icell->shg(xx,yy,1,inode);
								vars.at(variable_name).at(icell->t(inode)) += N * input_var.at(ii);
							break;
							case 4: //g2p
								N = icell->shp(xx, yy, inode);
								dprops.at(variable_name).at(ii) += N * input_var.at(icell->t(inode));
							break;
							case 5: //g2p dx
								N = icell->shg(xx,yy,0,inode);
								dprops.at(variable_name).at(ii) += N * input_var.at(icell->t(inode));
							break;
							case 6: //g2p dy
								N = icell->shg(xx,yy,1,inode);
								dprops.at(variable_name).at(ii) += N * input_var.at(icell->t(inode));
							break;
						}//switch
          }//inode for
        }//ii for
    }//icell for
    
    if (apply_mass)
    	switch ((transfer_code - 1) / 3){
    		case 0:
    			for (idx_t ii = 0; ii < M.size (); ++ii)
    				vars.at(variable_name).at(ii) /= M.at(ii);
    		break;
    		case 1:
    			for (idx_t ii = 0; ii < M.size (); ++ii)
    				dprops.at(variable_name).at(ii) /= M.at(ii);
    		break;
    	}
	}//transfer func

};//particles_t class

#endif /* PARTICLES_H */
