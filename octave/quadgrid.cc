// Copyright (C) 2022 Carlo de Falco
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program ; see the file COPYING.  If not, see
// <http://www.gnu.org/licenses/>.

#include "quadgrid.h"

DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA (quadgrid,
                                     "quadgrid",
                                     "quadgrid");

const std::shared_ptr<quadgrid_t<ColumnVector>>
ov_quadgrid (const octave_value &in) {
  const octave_base_value& rep = in.get_rep ();
  return static_cast<const quadgrid&> (rep).quadgrid_value ();
}

std::shared_ptr<quadgrid_t<ColumnVector>>
ov_quadgrid (octave_value &in) {
  const octave_base_value& rep = in.get_rep ();
  return static_cast<const quadgrid&> (rep).quadgrid_value ();
}

bool
is_quadgrid (const octave_value &in) {
  return (in.type_id () == quadgrid::static_type_id ());
}


DEFMETHOD_DLD (quadgrid, interp, args, ,"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} @var{QG} = quadgrid (@var{nx}, \
@var{hx}, @var{ny}, @var{hy})\n\
Return in @var{QG} a quadgrid object.\n\
@end deftypefn")
{
  octave_value retval;
  if (args.length () != 4)
    print_usage ();
  else
    {      
      retval = new quadgrid (args(0).idx_type_value (),
                             args(1).double_value (),
                             args(2).idx_type_value (),
                             args(3).double_value ());
    }  
  return retval;
}

// PKG_ADD: autoload ("quadgrid_loop", "quadgrid.oct");
DEFMETHOD_DLD (quadgrid_loop, interp, args, ,"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} @var{QG} = MPI_COMM_WORLD (@var{qg}, @var{nx}, \
@var{hx}, @var{ny}, @var{hy})\n\
Return in @var{QG} a quadgrid object.\n\
@end deftypefn")
{
  octave_value retval;
  if (args.length () != 1)
    print_usage ();
  else
    {    
      if (! is_quadgrid (args(0))) {
        error ("not a quadgrid");
      } else {
        auto grid = ov_quadgrid (args (0));
        std::cout << "rank " << grid->rank
                  << " owns " << grid->num_local_cells () << " cells"
                  << " of " << grid->num_global_cells () << " total,"
                  << " it owns " << grid->num_owned_nodes () << " nodes"
                  << " and touches " << grid->num_local_nodes () << " nodes"
                  << " of " << grid->num_global_nodes () << " total"
                  << std::endl;

        for (auto icell = grid->begin_cell_sweep ();
             icell != grid->end_cell_sweep (); ++icell) {

          std::cout << "visiting cell " << icell->get_local_cell_idx ()
                    << " row = " << icell->row_idx ()
                    << " col = " << icell->col_idx ()
                    << std::endl;

          for (auto inode = 0; inode < quadgrid_t<ColumnVector>::cell_t::nodes_per_cell; ++inode) {
            std::cout << "\tnode " << inode << " global index " << icell->gt (inode)
                      << " (rank) local index " << icell->t (inode)
                      << " coordinates " << icell->p (0, inode) << ", "
                      << icell->p (1, inode) << std::endl;

          }

          for (auto iedge = 0; iedge < quadgrid_t<ColumnVector>::cell_t::edges_per_cell; ++iedge)
            if (icell->e (iedge) != quadgrid_t<ColumnVector>::cell_t::NOT_ON_BOUNDARY)
              std::cout << "\tedge " << iedge << " of this cell is on boundary face "
                        << icell->e (iedge) << std::endl;
        }
      }
    }
  return retval;
}
