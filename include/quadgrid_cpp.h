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
#include "tbasisfun.h"



/// @brief 
/// @tparam distributed_vector 
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

  MPI_Comm          comm;
  int               rank;
  int               size;

// Degree
  static constexpr idx_t px=3;
  static constexpr idx_t py=3;
  // Regularity
  static constexpr idx_t rx=2;
  static constexpr idx_t ry=2;


private:

  grid_properties_t grid_properties;

  mutable cell_t   current_cell;
  mutable cell_t   current_neighbor;
  

public:

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
/*
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
  }*/

/*
// Return global index of the BSpline coefficient 
// This version takes span indexes
  static idx_t gt (idx_t inode, idx_t row_span_idx, idx_t col_span_idx, idx_t num_rows) {
    if (inode <0 || inode >= (px+1)*(py+1))
    return -1;// or std::assert

    idx_t r1=inode%(py+1);
    idx_t c1=inode/(py+1);
    return global_BS_idx(row_span_idx-py+r1,col_span_idx-px+c1,num_rows);

  }
*/

// Return global index of the BSpline coefficient 
//This version takes the row and col indexes of the cell
  static idx_t gt (idx_t inode, idx_t ridx, idx_t cidx, idx_t num_rows,idx_t num_cols) const {
    if (inode <0 || inode >= (px+1)*(py+1))
    return -1;// or std::assert

    idx_t ii=cell2span(ridx,py,ry,num_rows);
    idx_t jj=cell2span(cidx,px,rx,num_cols);

    idx_t r1=inode%(py+1);
    idx_t c1=inode/(py+1);
    return global_BS_idx(ii-py+r1,jj-px+c1);

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

// Bsplines

// Restituisce l'indice della funzione di base dato l'indice di cella
// N è il numero di celle in quella direz=num breakpts-1 (numrows/numcols)
// p degree
// r regolarità

static idx_t cell2span( idx_t index, idx_t p, idx_t r, idx_t N){// index è indice di riga/col
std::assert(index>=0 && index<=N-1);

return p + (p-r)*index;
}

static idx_t knot2gind( idx_t k_ind, idx_t p, idx_t r, idx_t N){// N=nrows+1 o ncols+1
m=(N-2)*(p-r)+2*(p+1)-1;// max index of knot
std::assert(k_ind>=0 && k_ind<=m);

if(k_ind<=p){
  return 0;
}
if(k_ind>=m-p){//controllare coerenza con onebasisfun, m-p-1 ritorna N-2 senza pb
return N-1;
}
else{ 
  if ((k_ind-p)%(p-r)==0){
    return (k_ind-p)/(p-r);
  }
  else{ 
    return (k_ind-p)/(p-r)+1;
  } 
}
}

// Convert from 2d-basisfun indexes to global BSpline index
// NB Number of basis functions is less than number of knots

static idx_t global_BS_idx(idx_t ii, idx_t jj, idx_t num_rows){// num_rows è il numero di celle, a me serve il numero di breakpts-> +1
return ii+jj*((num_rows-1)*(px-rx)+(px+1)*2-px-1);
}

struct
  Span_iterator {
    using difference_type = typename std::make_signed_t<idx_t>;
    using iterator_category = std::random_access_iterator_tag;
    using value_type = idx_t;
    using reference = idx_t;
    using pointer = const idx_t*;
    iterator &operator ++() { ++i_; return *this; }
    iterator operator ++(int) { iterator copy(*this); ++i_; return copy; }
    iterator &operator --() { --i_; return *this; }
    iterator operator --(int) { iterator copy(*this); --i_; return copy; }
    iterator& operator +=(difference_type n) { i_ = static_cast<value_type>(static_cast<difference_type>(i_) + n); return *this; }
    iterator& operator -=(difference_type n) { i_ = static_cast<value_type>(static_cast<difference_type>(i_) - n); return *this; }


    double operator *() const {
      return quadgird_t::knot2gind(i_,px,rx,grid_properties.numcols+1)*grid_properties.hx;
 }
    bool operator ==(const iterator &other) const { return i_ == *other; }
    bool operator !=(const iterator &other) const { return i_ != *other; }

    bool operator <(const iterator &other) const { return i_ < *other; }
    difference_type operator -(const iterator &other) const { return static_cast<difference_type>(i_) - static_cast<difference_type>(*other); }
    iterator operator -(const difference_type other) const { return iterator (static_cast<value_type> (static_cast<difference_type>(i_) - other)); }
    iterator operator +(const difference_type other) const { return iterator (static_cast<value_type> (static_cast<difference_type>(i_) + other)); }
    value_type operator[] (const value_type& idx) { return i_ + idx; }
    protected: explicit iterator(difference_type start) : i_ (start) {};
  private:
    idx_t i_;
  };
  



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
  shp (double x, double y, idx_t inode,
       idx_t c, idx_t r, double hx, double hy) {

    if (inode <0 || inode >= (px+1)*(py+1))
      return -1;// or std::assert

    idx_t ii=cell2span(r,py,ry,num_rows);
    idx_t jj=cell2span(c,px,rx,num_cols);

    idx_t r1=inode%(py+1);
    idx_t c1=inode/(py+1);
    

    return onebasisfun2d (x, y, px, py, Ubegin, Uend, Vbegin, Vend);
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

    static constexpr idx_t nodes_per_cell = 16;
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
	    return quadgrid_t::gt (i, row_idx (), col_idx (), num_rows (), num_cols());
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



};



#include "quadgrid_cpp_imp.h"

#endif /* QUADGRID_H */

