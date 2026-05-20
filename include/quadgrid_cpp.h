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
    idx_t             px, py;  // Degree
    idx_t             rx, ry;  // Regularity
    idx_t             nodes_per_cell;  // (px+1)*(py+1)
    idx_t             num_dof_x;
    idx_t             num_dof_y;
    std::vector<double> knot_vect_x;
    std::vector<double> knot_vect_y;

    idx_t             start_cell_row;
    idx_t             end_cell_row;
    idx_t             start_cell_col;
    idx_t             end_cell_col;
    idx_t             start_owned_nodes;
    idx_t             num_owned_nodes;
    
  };

  std::vector<double> Init_knot_vector(double h, idx_t n, idx_t p, idx_t r) const{
  std::vector<double> grid_nodes(n+1);

  for(idx_t i=0; i<=n;++i)
  grid_nodes[i]=i*h;

  return bspline::open_knot_vector(grid_nodes.cbegin(),grid_nodes.cend(),p,r);
 }


  void
  from_json (const nlohmann::json &j, grid_properties_t &q) {

    j.at ("nx").get_to (q.numcols);
    j.at ("ny").get_to (q.numrows);
    j.at ("hx").get_to (q.hx);
    j.at ("hy").get_to (q.hy);
    q.px = j.value("px", 3);  
    q.py = j.value("py", 3);  
    q.rx = j.value("rx", 2);  
    q.ry = j.value("ry", 2);  

    q.nodes_per_cell = (q.px + 1) * (q.py + 1);
    q.num_dof_x = (q.numcols - 1) * (q.px - q.rx) + 2 * (q.px + 1) - q.px - 1;
    q.num_dof_y = (q.numrows - 1) * (q.py - q.ry) + 2 * (q.py + 1) - q.py - 1;
    q.start_cell_row = 0;
    q.end_cell_row = q.numrows - 1;
    q.start_cell_col = 0;
    q.end_cell_col = q.numcols - 1;
    q.start_owned_nodes = 0;
    q.num_owned_nodes = q.num_dof_y * q.num_dof_x; //forse meglio cambiare nodes con dof

    q.knot_vect_y = Init_knot_vector(q.hy,q.numrows,q.py,q.ry);
    q.knot_vect_x = Init_knot_vector(q.hx,q.numcols,q.px,q.rx);
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
  static idx_t gt (idx_t inode, idx_t col_span_idx, idx_t row_span_idx, idx_t num_rows) {
    if (inode <0 || inode >= (px+1)*(py+1))
    return -1;// or std::assert

    idx_t r1=inode%(py+1);
    idx_t c1=inode/(py+1);
    return sub2gind(row_span_idx-py+r1,col_span_idx-px+c1,((num_rows-1)*(py-ry)+(py+1)*2-py-1);

  }
*/
// Posso sostituirci N_dof_y migliorandola, ma viene usato dai g2p_helper_t che potrei rimuovere
// Return global index of the BSpline coefficient 
//This version takes the row and col indexes of the cell
  static idx_t gt (idx_t inode, idx_t cidx, idx_t ridx, idx_t N_dof_y, idx_t px, idx_t py, idx_t rx, idx_t ry) {
    if (inode <0 || inode >= (px+1)*(py+1))
    return -1;// or std::assert

    idx_t ii=cell2span(ridx,py,ry);
    idx_t jj=cell2span(cidx,px,rx);

    idx_t r1=inode%(py+1);
    idx_t c1=inode/(py+1);
    return sub2gind(ii-py+r1,jj-px+c1,N_dof_y);// Third param is N_dof_y 

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
/*
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
  }*/

// Bsplines

// Restituisce l'indice della funzione di base dato l'indice di cella (riga o col)
// p degree
// r regolarità
static idx_t cell2span( idx_t index, idx_t p, idx_t r){ 
// std::assert(index>=0 && index<=N-1);

return p + (p-r)*index;
}

/*
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
*/
// Convert from 2d-basisfun indexes to global BSpline index
// NB Number of basis functions is less than number of knots

// static idx_t global_Spline_idx(idx_t ii, idx_t jj, idx_t num_rows){// num_rows è il numero di celle, a me serve il numero di breakpts-> +1
// return ii+jj*((num_rows-1)*(py-ry)+(py+1)*2-py-1);
// }

/*
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
    
  }*/

 static double
  shp (double x, double y, idx_t inode,
       idx_t c, idx_t r, idx_t px, idx_t py, 
       idx_t rx, idx_t ry, idx_t num_dof_x, idx_t num_dof_y, std::vector<double> const  & knot_vect_x, std::vector<double> const & knot_vect_y) {
        using namespace bspline;
    if (inode <0 || inode >= (px+1)*(py+1))
      return -1;// or std::assert

    idx_t ii=cell2span(r,py,ry);
    idx_t jj=cell2span(c,px,rx);

    idx_t r1=inode%(py+1);
    idx_t c1=inode/(py+1);
    

    std::vector<double>::const_iterator Span_x_begin = knot_vect_x.begin() + jj-px+c1;
    std::vector<double>::const_iterator Span_x_end = knot_vect_x.begin() + jj+c1+2;
    std::vector<double>::const_iterator Span_y_begin = knot_vect_y.begin() + ii-py+r1;
    std::vector<double>::const_iterator Span_y_end = knot_vect_y.begin() + ii+r1+2;

    if (ii-py+r1==num_dof_y-1){// check if the bspline index is the last index possible, i.e. we are on the boundary
      if (jj-px+c1==num_dof_x-1)
        return onebasisfun2d<Position::Boundary,Position::Boundary>
         (x, y, px, py, Span_x_begin, Span_x_end, Span_y_begin, Span_y_end);
      else
        return onebasisfun2d<Position::Internal,Position::Boundary> (x, y, px, py, Span_x_begin, Span_x_end, Span_y_begin, Span_y_end);
    }
    if (jj-px+c1==num_dof_x-1)
      return onebasisfun2d<Position::Boundary,Position::Internal> (x, y, px, py,
         Span_x_begin, Span_x_end, Span_y_begin, Span_y_end);
    return onebasisfun2d<Position::Internal,Position::Internal> (x, y, px, py,
       Span_x_begin, Span_x_end, Span_y_begin, Span_y_end);
       }


  static double
  shg (double x, double y, idx_t idir, idx_t inode,
       idx_t c, idx_t r, idx_t px, idx_t py, idx_t rx, idx_t ry,
      idx_t num_dof_x, idx_t num_dof_y, std::vector<double> const  & knot_vect_x, std::vector<double> const & knot_vect_y) {
   using namespace bspline;
    if (inode <0 || inode >= (px+1)*(py+1))
      return -1;// or std::assert
    //std::assert(idir==0 || idir==1)

    idx_t ii=cell2span(r,py,ry);
    idx_t jj=cell2span(c,px,rx);
        
    idx_t r1=inode%(py+1);
    idx_t c1=inode/(py+1);
    
  
    std::vector<double>::const_iterator Span_x_begin = knot_vect_x.begin() + jj-px+c1;
    std::vector<double>::const_iterator Span_x_end = knot_vect_x.begin() + jj+c1+2;
    std::vector<double>::const_iterator Span_y_begin = knot_vect_y.begin() + ii-py+r1;
    std::vector<double>::const_iterator Span_y_end = knot_vect_y.begin() + ii+r1+2;


if (idir==0){//x-deriv
if (ii-py+r1==num_dof_y-1){// check if the bspline index is the last index possible, i.e. we are on the boundary
      if (jj-px+c1==num_dof_x-1)
        return onebasisfun<Position::Boundary>(y,py,Span_y_begin, Span_y_end)*onebasisfunder<Position::Boundary>(x,px,Span_x_begin, Span_x_end);
      else
        return onebasisfun<Position::Boundary>(y,py,Span_y_begin, Span_y_end)*onebasisfunder<Position::Internal>(x,px,Span_x_begin, Span_x_end);
    }
if (jj-px+c1==num_dof_x-1)
      return onebasisfun<Position::Internal>(y,py,Span_y_begin, Span_y_end)*onebasisfunder<Position::Boundary>(x,px,Span_x_begin, Span_x_end);
return onebasisfun<Position::Internal>(y,py,Span_y_begin, Span_y_end)*onebasisfunder<Position::Internal>(x,px,Span_x_begin, Span_x_end);
}
else if(idir==1){//y-deriv
  if (ii-py+r1==num_dof_y-1){// check if the bspline index is the last index possible, i.e. we are on the boundary
      if (jj-px+c1==num_dof_x-1)
        return onebasisfunder<Position::Boundary>(y,py,Span_y_begin, Span_y_end)*onebasisfun<Position::Boundary>(x,px,Span_x_begin, Span_x_end);
      else
        return onebasisfunder<Position::Boundary>(y,py,Span_y_begin, Span_y_end)*onebasisfun<Position::Internal>(x,px,Span_x_begin, Span_x_end);
    }
if (jj-px+c1==num_dof_x-1)
      return onebasisfunder<Position::Internal>(y,py,Span_y_begin, Span_y_end)*onebasisfun<Position::Boundary>(x,px,Span_x_begin, Span_x_end);

return onebasisfunder<Position::Internal>(y,py,Span_y_begin, Span_y_end)*onebasisfun<Position::Internal>(x,px,Span_x_begin, Span_x_end);


}
return -1;
}



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


  protected :
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

  // private :
  //   cell_t *data;// Ci sarebbero 2 data così
  };




  class
  cell_t
  {

    friend class cell_iterator;
    friend class quadgrid_t;

  public:

    static constexpr idx_t edges_per_cell = 4;
    static constexpr idx_t NOT_ON_BOUNDARY = -1;

    cell_t (const grid_properties_t& _gp)
      : is_ghost (false), rowidx (0), colidx (0), grid_properties (_gp) { };

    //double
    //p (idx_t i, idx_t j) const;

    // double
    // centroid (idx_t i);

    // idx_t
    // t (idx_t i) const;

   idx_t
    gt (idx_t i) const {  
	    return quadgrid_t::gt (i, col_idx (), row_idx (), grid_properties.num_dof_y , grid_properties.px,
       grid_properties.py, grid_properties.rx, grid_properties.ry);
    }
  
    idx_t
    e (idx_t i) const;

    double
    shp (double x, double y, idx_t inode) const;

    // double
    // shp_new (double x, double y, idx_t inode) const;

    double
    shg (double x, double y, idx_t idir, idx_t inode) const;

    // neighbor_iterator
    // begin_neighbor_sweep ();

    // const neighbor_iterator
    // begin_neighbor_sweep () const;

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
    nodes_per_cell () const
    { return grid_properties.nodes_per_cell; };

    idx_t
    row_idx () const
    { return rowidx; };

    idx_t
    col_idx () const
    { return colidx; };
    
    // Potrebbero essere fuorvianti, capire dove serve indice totale cella
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
    grid_properties.px = 3;  
    grid_properties.py = 3;  
    grid_properties.rx = 2; 
    grid_properties.ry = 2;  
    grid_properties.num_dof_x = 0;
    grid_properties.num_dof_y = 0;
    grid_properties.nodes_per_cell = (grid_properties.px + 1) * (grid_properties.py + 1); 
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
             double hx, double hy);// Poi devo riaggiornare il knot_vector (fatto)

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
  num_owned_nodes () const
  { return grid_properties.num_owned_nodes; };

  idx_t
  num_local_nodes () const;

  idx_t
  num_global_nodes () const;

  idx_t
  num_local_cells () const;

  idx_t
  num_global_cells () const;

  idx_t N_dof_x() const
  {return grid_properties.num_dof_x;};

  idx_t N_dof_y() const
  {return grid_properties.num_dof_y;};


  std::vector<double> const & knot_vector_x() const {
    return grid_properties.knot_vect_x;
  };

  std::vector<double> const & knot_vector_y() const {
    return grid_properties.knot_vect_y;
  };


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
  nodes_per_cell () const
  { return grid_properties.nodes_per_cell; };

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

private:

  grid_properties_t grid_properties;

  mutable cell_t   current_cell;
  mutable cell_t   current_neighbor;

};



#include "quadgrid_cpp_imp.h"

#endif /* QUADGRID_H */

