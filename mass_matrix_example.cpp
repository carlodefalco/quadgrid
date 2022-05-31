#include <quadgrid_cpp.h>


int
main (int argc, char *argv[]) {

  MPI_Init (&argc, &argv);
  
  quadgrid_t<std::vector<double>> grid;
  grid.set_sizes (1000, 1000, .1, .2);
  std::vector<double> mass (grid.num_global_nodes (), 0.0);
  
  if (grid.rank == 0) {
    std::cout << "grid created with " << grid.num_rows ()
              << " rows and " << grid.num_cols () << " columns"
              << std::endl;
 
    for (auto icell = grid.begin_cell_sweep ();
         icell != grid.end_cell_sweep (); ++icell) {
      
      for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode) {
        mass[icell->gt (inode)] += (grid.hx () / 2.) * (grid.hy () / 2.);          
      }
    }     
  
    // uncomment to print all 1e6 elements of "mass"
    // for (auto ii : mass) {
    //   std::cout << ii << std::endl;
    // }

  }
  MPI_Finalize ();
  return 0;
};


