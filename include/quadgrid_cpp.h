#ifndef QUADGRID_H
#define QUADGRID_H

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <json.hpp>
#include <map>
#ifdef USE_MPI_H
#include <mpi.h>
#else
#define MPI_Comm int
#define MPI_COMM_WORLD 1
#define MPI_Initialized(x)
#define MPI_Comm_size(x, y) { *y = 0; }
#define MPI_Comm_rank(x, y) { *y = 0; }
#endif
#include <vector>


template <class distributed_vector>
class
quadgrid_t
{

public:

  using idx_t = int;

  class  cell_t;

  struct grid_properties_t {
    idx_t             numrows;
    idx_t             numcols;
    double            hx, hy;
    idx_t             start_cell_row;
    idx_t             end_cell_row;
    idx_t             start_cell_col;
    idx_t             end_cell_col;
    idx_t             start_owned_nodes;
    idx_t             num_owned_nodes;
  };

  void
  from_json (const nlohmann::json &j, grid_properties_t &q) {

    j.at ("nx").get_to (q.numcols);
    j.at ("ny").get_to (q.numrows);
    j.at ("hx").get_to (q.hx);
    j.at ("hy").get_to (q.hy);


    q.start_cell_row = 0;
    q.end_cell_row = q.numrows - 1;
    q.start_cell_col = 0;
    q.end_cell_col = q.numcols - 1;
    q.start_owned_nodes = 0;
    q.num_owned_nodes = (q.numrows+1)*(q.numcols+1);

  }

  static idx_t
  gind2col (idx_t idx, idx_t numrows) {
    return (idx / numrows);
  }

  static idx_t
  gind2row (idx_t idx, idx_t numrows) {
    return  (idx % numrows);
  }

  static idx_t
  gt (idx_t inode, idx_t cidx, idx_t ridx, idx_t numrows) {
    idx_t bottom_left = 0;
    // should check that inode < 4 in an efficient way
    bottom_left =  ridx + cidx * (numrows + 1);
    switch (inode) {
    case 0 :
      return (bottom_left);
      break;
    case 1 :
      return (bottom_left + 1);
      break;
    case 2 :
      return (bottom_left + (numrows + 1));
      break;
    case 3 :
      return (bottom_left + (numrows + 2));
      break;
    default :
      return -1;
    }
  }


  //-----------------------------------
  //
  //   Numbering of nodes and edges :
  //
  //              1
  //              |
  //              V
  // 1 -> O---------------O <- 3
  //      |               |
  //      |               |
  // 2 -> |               | <- 3
  //      |               |
  //      |               |
  // 0 -> O---------------O <- 2
  //              ^
  //              |
  //              0
  //
  //-----------------------------------

  static double
  p (idx_t idir, idx_t inode, idx_t colidx, idx_t rowidx, double hx, double hy) {
    double bottom_left = 0.0;
    // should check that inode < 4 in an efficient way
    if (idir == 0) {
      bottom_left = colidx * hx;
      if (inode > 1)
	bottom_left += hx;
    } else {
      bottom_left = rowidx * hy;
      if (inode == 1 || inode == 3)
	bottom_left += hy;
    }
    return (bottom_left);
  }

  
  static double
  shp (double x, double y, idx_t inode,
       idx_t c, idx_t r, double hx, double hy) {

    switch (inode) {
    case 3 :
      return ((x - p(0,0,c,r,hx,hy))/hx * (y - p(1,0,c,r,hx,hy))/hy);
      break;
    case 2 :
      return ((x - p(0,0,c,r,hx,hy))/hx * (1. - (y - p(1,0,c,r,hx,hy))/hy));
      break;
    case 1 :
      return ((1. - (x - p(0,0,c,r,hx,hy))/hx) * (y - p(1,0,c,r,hx,hy))/hy);
      break;
    case 0 :
      return ((1. - (x - p(0,0,c,r,hx,hy))/hx) * (1. - (y - p(1,0,c,r,hx,hy))/hy));
      break;
    default :
      return 0;
    }
    
  }

  static double
  shg (double x, double y, idx_t idir, idx_t inode,
       idx_t c, idx_t r, double hx, double hy) {
    switch (inode) {
    case 3 :
      if (idir == 0) {
	return ((1. / hx) * ((y - p(1,0,c,r,hx,hy)) / hy));
      }
      else if (idir == 1) {
	return (((x - p(0,0,c,r,hx,hy)) / hx) * (1. / hy));
      }
      break;
    case 2 :
      if (idir == 0) {
	return ((1. / hx) * ((1. - (y - p(1,0,c,r,hx,hy)) / hy)));
      }
      else if (idir == 1) {
	return (((x - p(0,0,c,r,hx,hy)) / hx) * (- 1. / hy));
      }
      break;
    case 1 :
      if (idir == 0) {
	return ((- 1. / hx) * ((y - p(1,0,c,r,hx,hy)) / hy));
      }
      else if (idir == 1) {
	return ((1. - (x - p(0,0,c,r,hx,hy)) / hx) * (1. / hy));
      }
      break;
    case 0 :
      if (idir == 0) {
	return ((- 1. / hx) * (1. - (y - p(1,0,c,r,hx,hy)) / hy));
      }
      else if (idir == 1) {
	return ((1. - (x - p(0,0,c,r,hx,hy))/ hx) * (- 1. / hy));
      }
      break;
    }
    return 0.;
  };

  
  static idx_t
  sub2gind (idx_t r, idx_t c, idx_t nr) {
    return  (r + nr * c);
  }
  
  
  class
  cell_iterator
  {

  public:

    cell_iterator (cell_t* _data = nullptr)
      : data (_data) { };

    void
    operator++ ();

    cell_t&
    operator* ()
    { return *(this->data); };

    const cell_t&
    operator* () const
    { return *(this->data); };

    cell_t*
    operator-> ()
    { return this->data; };

    const cell_t*
    operator-> () const
    { return this->data; };

    bool
    operator== (const cell_iterator& other)
    { return (data == other.data); }

    bool
    operator!= (const cell_iterator& other)
    { return ! ((*this) == other); }


  private :
    cell_t *data;
  };

  class
  neighbor_iterator : public cell_iterator
  {

  public:

    void
    operator++ ();

    neighbor_iterator (cell_t *_data = nullptr,
                       int _face_idx = -1)
      : cell_iterator (_data), face_idx (_face_idx) { };

    int
    get_face_idx ()
    { return face_idx; };

  private:
    idx_t face_idx; /// Face index in 0...3 (-1 if not defined).

  private :
    cell_t *data;
  };

  class
  cell_t
  {

    friend class cell_iterator;
    friend class quadgrid_t;

  public:

    static constexpr idx_t nodes_per_cell = 4;
    static constexpr idx_t edges_per_cell = 4;
    static constexpr idx_t NOT_ON_BOUNDARY = -1;

    cell_t (const grid_properties_t& _gp)
      : grid_properties (_gp), rowidx (0), colidx (0), is_ghost (false) { };

    double
    p (idx_t i, idx_t j) const;

    double
    centroid (idx_t i);

    idx_t
    t (idx_t i) const;

    idx_t
    gt (idx_t i) const {  
      // should check that inode < 4 in an efficient way
      switch (i) {
      case 0 :
      case 1 :
      case 2 :
      case 3 :
	return quadgrid_t::gt (i, col_idx (), row_idx (), num_rows ());
	break;
      default :
	return -1;
      }
    }
  
    idx_t
    e (idx_t i) const;

    double
    shp (double x, double y, idx_t inode) const;

    double
    shp_new (double x, double y, idx_t inode) const;

    double
    shg (double x, double y, idx_t idir, idx_t inode) const;

    neighbor_iterator
    begin_neighbor_sweep ();

    const neighbor_iterator
    begin_neighbor_sweep () const;

    neighbor_iterator
    end_neighbor_sweep ()
    { return neighbor_iterator (); };

    const neighbor_iterator
    end_neighbor_sweep () const
    { return neighbor_iterator (); };

    idx_t
    get_local_cell_idx () const
    { return local_cell_idx; };

    idx_t
    get_global_cell_idx () const
    { return global_cell_idx; };

    idx_t
    end_cell_col () const
    { return grid_properties.end_cell_col; };

    idx_t
    end_cell_row () const
    { return grid_properties.end_cell_row; };

    idx_t
    start_cell_col () const
    { return grid_properties.start_cell_col; };

    idx_t
    start_cell_row () const
    { return grid_properties.start_cell_row; };

    idx_t
    num_rows () const
    { return grid_properties.numrows; };

    idx_t
    num_cols () const
    { return grid_properties.numcols; };

    idx_t
    row_idx () const
    { return rowidx; };

    idx_t
    col_idx () const
    { return colidx; };
    
    
    idx_t
    sub2gind (idx_t r, idx_t c) const {
      return  quadgrid_t::sub2gind (r, c, grid_properties.numrows);
    }

    idx_t
    gind2row (idx_t idx) const {
      return quadgrid_t:: gind2row (idx, grid_properties.numrows);
    }

    idx_t
    gind2col (idx_t idx) const {
      return quadgrid_t:: gind2col (idx, grid_properties.numrows);
    }

    void
    reset () {
      rowidx = grid_properties.start_cell_row;
      colidx = grid_properties.start_cell_col;
      global_cell_idx = sub2gind (rowidx, colidx);
      local_cell_idx = global_cell_idx -
        sub2gind (grid_properties.start_cell_row,
                  grid_properties.start_cell_col);
    };

  private:

    bool                     is_ghost;
    idx_t                    rowidx;
    idx_t                    colidx;
    idx_t                    local_cell_idx;
    idx_t                    global_cell_idx;
    const grid_properties_t &grid_properties;

  };


  /// Default constructor, set all pointers to nullptr.
  quadgrid_t (MPI_Comm _comm = MPI_COMM_WORLD) :
    comm (_comm), rank (0), size (1),
    current_cell (grid_properties),
    current_neighbor (grid_properties)
  {
    int flag = 0;
    MPI_Initialized (&flag);
    if (flag) {
      MPI_Comm_rank (comm, &rank);
      MPI_Comm_size (comm, &size);
    } else {
      rank = 0;
      size = 1;
    }
    grid_properties.numrows = 0;
    grid_properties.numcols = 0;
    grid_properties.hx = 0.;
    grid_properties.hy = 0.;
    grid_properties.start_cell_row = 0;
    grid_properties.end_cell_row = 0;
    grid_properties.start_cell_col = 0;
    grid_properties.end_cell_col = 0;
    grid_properties.start_owned_nodes = 0;
    grid_properties.num_owned_nodes = 0;
  };

  /// Ctor that reads grid properties from a json object.
  quadgrid_t (const nlohmann::json &j, MPI_Comm _comm = MPI_COMM_WORLD) :
    quadgrid_t(_comm) { from_json (j, grid_properties); };

  /// Delete copy constructor.
  quadgrid_t (const quadgrid_t &) = delete;

  /// Delete assignment operator.
  quadgrid_t &
  operator= (const quadgrid_t &) = delete;

  /// Destructor.
  ~quadgrid_t () = default;

  void
  set_sizes (idx_t numrows, idx_t numcols,
             double hx, double hy);

  void
  vtk_export (const char *filename,
              const std::map<std::string,
              distributed_vector> & f) const;

  void
  octave_ascii_export (const char *filename,
                       const std::map<std::string,
                       distributed_vector> & f) const;

  cell_iterator
  begin_cell_sweep ();

  const cell_iterator
  begin_cell_sweep () const;

  cell_iterator
  end_cell_sweep ()
  { return cell_iterator (); };

  const cell_iterator
  end_cell_sweep () const
  { return cell_iterator (); };

  idx_t
  num_owned_nodes ()
  { return grid_properties.num_owned_nodes; };

  idx_t
  num_local_nodes () const;

  idx_t
  num_global_nodes () const;

  idx_t
  num_local_cells () const;

  idx_t
  num_global_cells () const;

  idx_t
  num_rows () const
  { return grid_properties.numrows; };

  idx_t
  num_cols () const
  { return grid_properties.numcols; };

  double
  hx () const
  { return grid_properties.hx; };

  double
  hy () const
  { return grid_properties.hy; };

  idx_t
  sub2gind (idx_t r, idx_t c) const {
    return  (r + grid_properties.numrows * c);
  }

  idx_t
  gind2row (idx_t idx) const {
    return quadgrid_t:: gind2row (idx, grid_properties.numrows);
  }

  idx_t
  gind2col (idx_t idx) const {
    return quadgrid_t:: gind2col (idx, grid_properties.numrows);
  }

  const cell_t&
  operator[] (idx_t tmp) const;

  MPI_Comm          comm;
  int               rank;
  int               size;

private :

  grid_properties_t grid_properties;

  mutable cell_t   current_cell;
  mutable cell_t   current_neighbor;


};



#include "quadgrid_cpp_imp.h"

#endif /* QUADGRID_H */

