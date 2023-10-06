#ifndef HAVE_CUSTOM_ITERATOR_GRID_H
#define HAVE_CUSTOM_ITERATOR_GRID_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <json.hpp>
#include <map>
#include <vector>

constexpr int NOT_ON_BOUNDARY = -1;

struct grid_t;

struct grid_t {
  
  const int num_rows; /*!< number of cell rows */
  const int num_cols; /*!< number of cell columns */
  const double hx;    /*!< cell width */
  const double hy;    /*!< cell height */

  /**
   *  \brief cell class
   */
  struct cell_t {
    int row, column;
    const grid_t *grid; /*!< Pointer to the enclosing grid object */
  
    cell_t () = delete;
    cell_t (const grid_t *g_)
      : grid (g_) {};

    int
    index () const {
      return row + grid->num_rows * column;
    }

    void
    set_index (int idx) {
      column = idx / grid->num_rows;
      row = idx % grid->num_rows;
      if (column >= grid->num_cols) {
	column = -1;
	row = -1;
      }
    }

    int
    t (int inode) const {
      int r = row, c = column;
      if (inode < 0 || inode > 3) throw;
      if (inode == 1 || inode == 3) ++r;
      if (inode == 2 || inode == 3) ++c;
      return r + (grid->num_rows + 1) * c;
    }

    int
    gt (int inode) const {
      return t (inode);
    }
    
    double
    p (int dir, int inode) const {
      switch (dir) {
      case 0:
	if (inode == 0 || inode == 1)
	  return column * grid->hx;
	else if (inode == 2 || inode == 3)
	  return (column + 1) * grid->hx;
	break;
      case 1:
	if (inode == 0 || inode == 2)
	  return row * grid->hy;
	else if (inode == 1 || inode == 3)
	  return (row + 1) * grid->hy;
	break;
      default:
	throw;
      }
      return 0;
    }

    double
    shp (double x, double y, int inode) const {
      switch (inode) {
      case 3 :
	return ((x - p(0,0))/grid->hx *
		(y - p(1,0))/grid->hy);
	break;
      case 2 :
	return ((x - p(0,0))/grid->hx *
		(1. - (y - p(1,0))/grid->hy));
	break;
      case 1 :
	return ((1. - (x - p(0,0))/grid->hx) *
		(y - p(1,0))/grid->hy);
	break;
      case 0 :
	return ((1. - (x - p(0,0))/grid->hx) *
		(1. - (y - p(1,0))/grid->hy));
	break;
      default :
	throw std::out_of_range ("inode must be in range 0..3");
      }
    }
    
    double
    shg (double x, double y, int idir, int inode) const {
      switch (inode) {
      case 3 :
	if (idir == 0) {
	  return ((1. / grid->hx) *
		  ((y - p(1,0)) / grid->hy));
	}
	else if (idir == 1) {
	  return (((x - p(0,0)) / grid->hx) *
		  (1. / grid->hy));
	}
	break;
      case 2 :
	if (idir == 0) {
	  return ((1. / grid->hx) *
		  ((1. - (y - p(1,0)) / grid->hy)));
	}
	else if (idir == 1) {
	  return (((x - p(0,0)) / grid->hx) *
		  (- 1. / grid->hy));
	}
	break;
      case 1 :
	if (idir == 0) {
	  return ((- 1. / grid->hx) *
		  ((y - p(1,0)) / grid->hy));
	}
	else if (idir == 1) {
	  return ((1. - (x - p(0,0)) / grid->hx) *
		  (1. / grid->hy));
	}
	break;
      case 0 :
	if (idir == 0) {
	  return ((- 1. /grid->hx) *
		  (1. - (y - p(1,0)) / grid->hy));
	}
	else if (idir == 1) {
	  return ((1. - (x - p(0,0))/grid->hx) *
		  (- 1. / grid->hy));
	}
	break;
      default :
	throw std::out_of_range ("inode must be in range 0..3, idir must be either 0 or 1");
      }
      return 0.;
    }
       
    int
    e (int iedge) const {

      if (row == 0 && iedge == 0)
	return 0;

      if (row == grid->num_rows - 1 && iedge == 1)
	return 1;

      if (column == 0 && iedge == 2)
	return 2;

      if (column == grid->num_cols - 1 && iedge == 3)
	return 3;

      return (NOT_ON_BOUNDARY);
    }
    
    void
    print () const {
      std::cout << "cell number " << index ()
		<< ": row " << row << " of " << grid->num_rows
		<< ", column " << column << " of " << grid->num_cols
		<< ", vertices: "
		<< gt(0) << "(" << p(0, 0) << ", " << p(1, 0) << "), "
		<< gt(1) << "(" << p(0, 1) << ", " << p(1, 1) << "), "
		<< gt(2) << "(" << p(0, 2) << ", " << p(1, 2) << "), "
		<< gt(3) << "(" << p(0, 3) << ", " << p(1, 3) << ") "
		<< std::endl;
    }
    
  };
  
  struct iterator
  {
    
    using value_type = cell_t;
    using difference_type = int;
    using pointer = cell_t*;
    using reference = cell_t&;
    using iterator_category = std::input_iterator_tag;
    
    mutable cell_t buffer;
    grid_t *grid;
    
    iterator () = delete;
    iterator (grid_t *grid_)
      : grid(grid_), buffer(grid_) { }

    iterator&
    operator++ () {
      buffer.set_index (buffer.index () + 1);
      return *this;
    }

    bool
    operator==(const iterator& rhs) {
      bool test = false;
      test = (rhs.buffer.index () == buffer.index ());
      return test;
    }
    
    bool
    operator!=(const iterator& rhs) {
      return !(*this == rhs);
    }

    reference
    operator*() const {
      return buffer;
    }

  };
  
  grid_t (int num_rows_, int num_cols_, double hx_, double hy_)
    : num_rows(num_rows_), num_cols(num_cols_), hx(hx_), hy(hy_) {};

  
  iterator
  begin () {
    iterator it (this);
    it.buffer.row = 0;
    it.buffer.column = 0;
    return it;
  }

  iterator
  end () {
    iterator it (this);
    it.buffer.row = -1;
    it.buffer.column = -1;
    return it;
  }

  const iterator
  cbegin () const {
    return this->begin ();
  }

  const iterator
  cend () {
    return this->end ();
  }
  
  void
  vtk_export (const char *filename,
	      const std::map<std::string, std::vector<double>> & f) const {

    std::ofstream ofs (filename, std::ofstream::out);

    // This is the XML format of a VTS file to write :

    ofs <<
      "<VTKFile type=\"StructuredGrid\" version=\"StructuredGrid\" byte_order=\"LittleEndian\">\n\
 <StructuredGrid WholeExtent=\"0 " << num_rows << " 0 " << num_cols << " 0 0\">\n \
 <Piece Extent=\"0 " << num_rows << " 0 " << num_cols << " 0 0\">\n";

    ofs << "      <PointData Scalars=\"";
    for (auto const & ii : f) {
      ofs << ii.first << ",";
    }
    ofs  << "\">\n";

    for (auto const & ii : f) {
      ofs << "        <DataArray type=\"Float64\" Name=\"" << ii.first <<"\" format=\"ascii\">\n        ";
      for (auto const & jj : ii.second) {
	ofs << jj << " ";
      }
      ofs << std::endl << "        </DataArray>" << std::endl;
    }

    ofs << "      </PointData>\n";

    ofs << "      <Points>\n        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int ii = 0; ii <= num_cols; ++ii) {
      ofs << "          ";
      for (int jj = 0; jj <= num_rows; ++jj) {
	ofs << std::setprecision(16) << hx*ii << " " << hy*jj << " 0 ";
      }
      ofs << std::endl;
    }
    ofs << "        </DataArray>\n      </Points>\n";


    ofs <<
      "    </Piece>\n\
</StructuredGrid>\n\
</VTKFile>\n";

    ofs.close ();
  }

  void
  octave_ascii_export
  (const char *filename,
   std::map<std::string, std::vector<double>> const & vars) const {

    std::ofstream os (filename, std::ofstream::out);
  
    os << "# name: p" << std::endl
       << "# type: matrix" << std::endl
       << "# rows: 2" << std::endl
       << "# columns: " << (num_cols + 1)*(num_rows + 1)
       << std::endl;

    for (int jj = 0; jj < num_cols + 1; ++jj) {
      for (int ii = 0; ii < num_rows + 1; ++ii) {
	os  << std::setprecision(16) << jj*hx << " ";
      }
    }
    os << std::endl;

    for (int jj = 0; jj < num_cols + 1; ++jj) {
      for (int ii = 0; ii < num_rows + 1; ++ii) {
	os  << std::setprecision(16) << ii*hy << " ";
      }
    }
    os << std::endl;
  
    os << "# name: t" << std::endl
       << "# type: matrix" << std::endl
       << "# rows: 4" << std::endl
       << "# columns: " << num_cols*num_rows << std::endl;
    for (auto inode = 0;
	 inode < 4;
	 ++inode) {
      for (auto icell = this->begin ();
	   icell != this->end (); ++icell) {
	os  << std::setprecision(16) << icell->gt(inode) << " ";
      }
      os << std::endl;
    }
  

    os << "# name: vars" << std::endl
       << "# type: scalar struct" << std::endl
       << "# ndims: 2" << std::endl
       << "1 1" << std::endl
       << "# length: " << vars.size () << std::endl;


    for (auto const & ii : vars) {
      os << "# name: " << ii.first << std::endl
	 << "# type: matrix" << std::endl
	 << "# rows: 1" << std::endl
	 << "# columns: " << ii.second.size () << std::endl;
      for (auto const & kk : ii.second) {
	os << std::setprecision(16) << kk << " ";
      }
      os << std::endl;
    }
    os << std::endl;

    os.close ();
  }
  
};


grid_t *
grid_from_json (const nlohmann::json &j) {

  int numcols, numrows;
  double hx, hy;
  j.at ("nx").get_to (numcols);
  j.at ("ny").get_to (numrows);
  j.at ("hx").get_to (hx);
  j.at ("hy").get_to (hy);
  return new grid_t (numrows, numcols, hx, hy);
};


/*
int
main () {
  grid_t grid (5, 4, .1, .2);

  std::cout << "iterate with begin/end \n\n";
  for (auto i = grid.begin (); i != grid.end (); ++i) {
    (*i).print ();
  }

  std::cout << "\n\niterate with range \n\n";
  for (auto const &i : grid) {
    i.print ();
  }

  std::vector<int> v(grid.num_rows*grid.num_cols, 0);
  std::transform (grid.begin (), grid.end (), v.begin (), [] (auto ii) { return ii.index (); });

  std::cout << "\n\nvector of indices \n\n";
  for (auto const &i : v) {
    std::cout << i << std::endl;
  }
  
  return 0;
}
*/

#endif
