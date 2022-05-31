#include <quadgrid_cpp.h>


int
main (int argc, char *argv[]) {

  MPI_Init (&argc, &argv);
  
  quadgrid_t<std::vector<double>> grid;
  grid.set_sizes (10, 5, .1, .2);

  if (grid.rank == 0) {
    std::cout << "grid created with " << grid.num_rows ()
              << " rows and " << grid.num_cols () << " columns"
              << std::endl;
  };

  for (auto irank = 0; irank < grid.size; ++irank) {
    if (irank == grid.rank) {
      std::cout << "rank " << grid.rank
                << " owns " << grid.num_local_cells () << " cells"
                << " of " << grid.num_global_cells () << " total,"
                << " it owns " << grid.num_owned_nodes () << " nodes"
                << " and touches " << grid.num_local_nodes () << " nodes"
                << " of " << grid.num_global_nodes () << " total"
                << std::endl;

      for (auto icell = grid.begin_cell_sweep ();
           icell != grid.end_cell_sweep (); ++icell) {

        std::cout << "visiting cell " << icell->get_local_cell_idx ()
                  << " row = " << icell->row_idx ()
                  << " col = " << icell->col_idx ()
                  << std::endl;

        for (auto inode = 0; inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell; ++inode) {
          std::cout << "\tnode " << inode << " global index " << icell->gt (inode)
                    << " (rank) local index " << icell->t (inode)
                    << " coordinates " << icell->p (0, inode) << ", "
                    << icell->p (1, inode) << std::endl;

        }

        for (auto iedge = 0; iedge < quadgrid_t<std::vector<double>>::cell_t::edges_per_cell; ++iedge)
          if (icell->e (iedge) != quadgrid_t<std::vector<double>>::cell_t::NOT_ON_BOUNDARY)
            std::cout << "\tedge " << iedge << " of this cell is on boundary face "
                      << icell->e (iedge) << std::endl;
        
          /*

          //
          // TO BE IMPLEMENTED : LOOP ON EDGE/CORNER NEIGHBOURS OF A CELL
          //

            for (auto jcell = icell->begin_neighbor_sweep ();
            icell != icell->end_neighbor_sweep (); ++icell) {

            }

          */
      }
      
    }
    MPI_Barrier (grid.comm);
  }
    

  MPI_Finalize ();
  return 0;
};


