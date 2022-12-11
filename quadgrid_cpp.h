#ifndef QUADGRID_H
#define QUADGRID_H

#include <fstream>
#include <mpi.h>
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
    gt (idx_t i) const;

    idx_t
    e (idx_t i) const;

    double
    shp (double x, double y, idx_t inode) const;

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
      return  (r + grid_properties.numrows * c);
    }

    idx_t
    gind2row (idx_t idx) const {
      return  (idx / grid_properties.numrows);
    }

    idx_t
    gind2col (idx_t idx) const {
      return  (idx % grid_properties.numrows);
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
              const distributed_vector & f) const;

  void
  vtk_export_cell (const char * filename,
                   const std::vector<double> & f) const;

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
    return  (idx / grid_properties.numrows);
  }

  idx_t
  gind2col (idx_t idx) const {
    return  (idx % grid_properties.numrows);
  }

  MPI_Comm          comm;
  int               rank;
  int               size;

private :

  mutable cell_t   current_cell;
  mutable cell_t   current_neighbor;

  grid_properties_t grid_properties;

};


template <class T>
typename quadgrid_t<T>::cell_iterator
quadgrid_t<T>::begin_cell_sweep () {
  current_cell.reset ();
  return cell_iterator (&current_cell);
}

template <class T>
const typename quadgrid_t<T>::cell_iterator
quadgrid_t<T>::begin_cell_sweep () const {
  current_cell.reset ();
  return cell_iterator (&current_cell);
}

template <class T>
void
quadgrid_t<T>::set_sizes (idx_t numrows, idx_t numcols,
                          double hx, double hy) {
  grid_properties.numrows = numrows;
  grid_properties.numcols = numcols;
  grid_properties.hx = hx;
  grid_properties.hy = hy;
  grid_properties.start_cell_row = 0;
  grid_properties.end_cell_row = numrows - 1;
  grid_properties.start_cell_col = 0;
  grid_properties.end_cell_col = numcols - 1;
  grid_properties.start_owned_nodes = 0;
  grid_properties.num_owned_nodes = (numrows+1)*(numcols+1);
}

template <class T>
void
quadgrid_t<T>::cell_iterator::operator++ () {
  static idx_t tmp;
  if (data != nullptr) {
    tmp =  data->rowidx + data->num_rows () *  data->colidx + 1;
    if (tmp > (data->end_cell_row () +
               data->num_rows () * data->end_cell_col ()))
      data = nullptr;
    else {
      data->rowidx = tmp % data->num_rows ();
      data->colidx = tmp / data->num_rows ();
      data->global_cell_idx = tmp;
      data->local_cell_idx = tmp -
        (data->start_cell_row () +
         data->num_rows () * data->start_cell_col ());
    }
  }
};

template <class T>
void
quadgrid_t<T>::neighbor_iterator::operator++ () {
  static idx_t tmp;
  if (data != nullptr) {
    tmp++;
    if (tmp >= quadgrid_t<T>::cell_t::edges_per_cell) {
      data = nullptr;
      face_idx = -1;
    }
    else {
      switch (face_idx) {
      case 0 :
        if (data->rowidx == 0)
          face_idx++;
        else {

          
            }
      case 1 :
        if (data->rowidx == data->num_rows () - 1)
          face_idx++;
        else {
        }
      case 2 :
        if (data->colidx == 1)
          face_idx++;
        else {
        }
      case 3 :
        if (data->colidx == data->num_cols () - 1) {
          face_idx = -1;
          data = nullptr;
        }
        else {
        }
      }
      data->rowidx = tmp % data->num_rows ();
      data->colidx = tmp / data->num_rows ();
      data->global_cell_idx = tmp;
      data->local_cell_idx = tmp -x
        (data->start_cell_row () +
         data->num_rows () * data->start_cell_col ());
    }
  }
};

template <class T>
typename quadgrid_t<T>::idx_t
quadgrid_t<T>::num_local_cells () const {
  return grid_properties.numrows*grid_properties.numcols;
}

template <class T>
typename quadgrid_t<T>::idx_t
quadgrid_t<T>::num_global_cells () const {
  return grid_properties.numrows*grid_properties.numcols;
}

template <class T>
typename quadgrid_t<T>::idx_t
quadgrid_t<T>::num_global_nodes () const {
  return (grid_properties.numrows+1)*(grid_properties.numcols+1);
}

template <class T>
typename quadgrid_t<T>::idx_t
quadgrid_t<T>::num_local_nodes () const {
  return (grid_properties.numrows+1)*(grid_properties.numcols+1);
}

template <class T>
typename quadgrid_t<T>::idx_t
quadgrid_t<T>::cell_t::gt (typename quadgrid_t<T>::idx_t inode) const {
  static idx_t bottom_left = 0;
  // should check that inode < 4 in an efficient way
  bottom_left =  row_idx () + col_idx () * (num_rows () + 1);
  switch (inode) {
  case 0 :
    return (bottom_left);
    break;
  case 1 :
    return (bottom_left + 1);
    break;
  case 2 :
    return (bottom_left + (num_rows () + 1));
    break;
  case 3 :
    return (bottom_left + (num_rows () + 2));
    break;
  default :
    return -1;
  }
}

template <class T>
typename quadgrid_t<T>::idx_t
quadgrid_t<T>::cell_t::t (typename quadgrid_t<T>::idx_t inode) const {
  static idx_t glob = 0;
  glob = gt (inode);
  // should check that inode < 4 in an efficient way
  if (glob < grid_properties.start_owned_nodes
      || glob >= (grid_properties.start_owned_nodes
                  + grid_properties.num_owned_nodes))
    return (glob);
  else
    return (glob - grid_properties.start_owned_nodes);
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


template <class T>
double
quadgrid_t<T>::cell_t::p (typename quadgrid_t<T>::idx_t idir,
                          typename quadgrid_t<T>::idx_t inode) const {
  static double bottom_left = 0.0;
  // should check that inode < 4 in an efficient way
  if (idir == 0) {
    bottom_left = col_idx () * grid_properties.hx;
    if (inode > 1)
      bottom_left += grid_properties.hx;
  } else {
    bottom_left = row_idx () * grid_properties.hy;
    if (inode == 1 || inode == 3)
      bottom_left += grid_properties.hy;
  }
  return (bottom_left);
}

template <class T>
typename quadgrid_t<T>::idx_t
quadgrid_t<T>::cell_t::e (typename quadgrid_t<T>::idx_t iedge) const {
  // should check that iedge < 4 in an efficient way

  if (row_idx () == 0 && iedge == 0)
    return 0;

  if (row_idx () == num_rows () - 1 && iedge == 1)
    return 1;

  if (col_idx () == 0 && iedge == 2)
    return 2;

  if (col_idx () == num_cols () - 1 && iedge == 3)
    return 3;

  return (NOT_ON_BOUNDARY);
}

template <class T>
double
quadgrid_t<T>::cell_t::shp (double x, double y, idx_t inode) const {
  switch (inode) {
  case 0 :
    return ((x - p(0,0))/grid_properties.hx *
            (y - p(1,0))/grid_properties.hy);
    break;
  case 1 :
    return ((x - p(0,0))/grid_properties.hx *
            (1. - (y - p(1,0))/grid_properties.hy));
    break;
  case 2 :
    return ((1. - (x - p(0,0))/grid_properties.hx) *
            (y - p(1,0))/grid_properties.hy);
    break;
  case 3 :
    return ((1. - (x - p(0,0))/grid_properties.hx) *
            (1. - (y - p(1,0))/grid_properties.hy));
    break;
  default :
    throw std::out_of_range ("inode must be in range 0..3");
  }
};

template <class T>
double
quadgrid_t<T>::cell_t::shg (double x, double y, idx_t idir, idx_t inode) const {
   switch (inode) {
  case 0 :
    if (idir == 0) {
      return ((1. / grid_properties.hx) *
              ((y - p(1,0)) / grid_properties.hy));
    }
    else if (idir == 1) {
      return (((x - p(0,0)) / grid_properties.hx) *
              (1. / grid_properties.hy));
    }    
  case 1 :
    if (idir == 0) {
      return ((1. / grid_properties.hx) *
              ((1. - (y - p(1,0)) / grid_properties.hy)));
    }
    else if (idir == 1) {
      return (((x - p(0,0)) / grid_properties.hx) *
              (- 1. / grid_properties.hy));
    }
  case 2 :
    if (idir == 0) {
    return ((- 1. / grid_properties.hx) *
            ((y - p(1,0)) / grid_properties.hy));
    }
    else if (idir == 1) {
      return ((1. - (x - p(0,0)) / grid_properties.hx) *
              (1. / grid_properties.hy));
    }
  case 3 :
    if (idir == 0) {
      return ((- 1. /grid_properties.hx) *
              (1. - (y - p(1,0)) / grid_properties.hy));
    }
    else if (idir == 1) {
    return ((1. - (x - p(0,0))/grid_properties.hx) *
            (- 1. / grid_properties.hy));
    }
  default :
    throw std::out_of_range ("inode must be in range 0..3, idir must be either 0 or 1");
  }
};



template <class T>
void
quadgrid_t<T>::vtk_export (const char *filename,
                           const T & f) const {

  std::ofstream ofs (filename, std::ofstream::out);

  // This is the XML format of a VTS file to write :

  /* <VTKFile type="StructuredGrid" ...> */
  /*   <StructuredGrid WholeExtent="x1 x2 y1 y2 z1 z2"> */
  /*   <Piece Extent="x1 x2 y1 y2 z1 z2"> */
  /*   <PointData>...</PointData> */
  /*   <CellData>...</CellData> */
  /*   <Points>...</Points> */
  /*   </Piece> */
  /*   </StructuredGrid> */
  /*   </VTKFile> */

};


#endif /* QUADGRID_H */

