
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
const typename quadgrid_t<T>::cell_t&
quadgrid_t<T>::operator[] (idx_t tmp) const {
  assert (tmp <= (current_cell.end_cell_row () +
                  current_cell.num_rows () * current_cell.end_cell_col ()));
  current_cell.rowidx = tmp % current_cell.num_rows ();
  current_cell.colidx = tmp / current_cell.num_rows ();
  current_cell.global_cell_idx = tmp;
  current_cell.local_cell_idx = tmp -
    (current_cell.start_cell_row () +
     current_cell.num_rows () * current_cell.start_cell_col ());
  return current_cell;
};



template <class T>
void
quadgrid_t<T>::neighbor_iterator::operator++ () {
  static idx_t tmp = 0; // Eseguita solo la prima volta
  if (data != nullptr) {
    tmp = ++face_idx;
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
	  data->rowidx = tmp % data->num_rows ();
	  data->colidx = tmp / data->num_rows ();
	  data->global_cell_idx = tmp;
	  data->local_cell_idx = tmp -
	    (data->start_cell_row () +
	     data->num_rows () * data->start_cell_col ());
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
      data->local_cell_idx = tmp -x // Probabile typo
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
  return grid_properties.numrows*grid_properties.numcols; // è uguale al locale
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





// Mi pare fosse già definita
/*template <class T>
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


*/


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
  case 3 :
    return ((x - p(0,0))/grid_properties.hx *
	    (y - p(1,0))/grid_properties.hy);
    break;
  case 2 :
    return ((x - p(0,0))/grid_properties.hx *
	    (1. - (y - p(1,0))/grid_properties.hy));
    break;
  case 1 :
    return ((1. - (x - p(0,0))/grid_properties.hx) *
	    (y - p(1,0))/grid_properties.hy);
    break;
  case 0 :
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
  case 3 :
    if (idir == 0) {
      return ((1. / grid_properties.hx) *
	      ((y - p(1,0)) / grid_properties.hy));
    }
    else if (idir == 1) {
      return (((x - p(0,0)) / grid_properties.hx) *
	      (1. / grid_properties.hy));
    }
    break;
  case 2 :
    if (idir == 0) {
      return ((1. / grid_properties.hx) *
	      ((1. - (y - p(1,0)) / grid_properties.hy)));
    }
    else if (idir == 1) {
      return (((x - p(0,0)) / grid_properties.hx) *
	      (- 1. / grid_properties.hy));
    }
    break;
  case 1 :
    if (idir == 0) {
      return ((- 1. / grid_properties.hx) *
	      ((y - p(1,0)) / grid_properties.hy));
    }
    else if (idir == 1) {
      return ((1. - (x - p(0,0)) / grid_properties.hx) *
	      (1. / grid_properties.hy));
    }
    break;
  case 0 :
    if (idir == 0) {
      return ((- 1. /grid_properties.hx) *
	      (1. - (y - p(1,0)) / grid_properties.hy));
    }
    else if (idir == 1) {
      return ((1. - (x - p(0,0))/grid_properties.hx) *
	      (- 1. / grid_properties.hy));
    }
    break;
  default :
    throw std::out_of_range ("inode must be in range 0..3, idir must be either 0 or 1");
  }
  return 0.;
};





template <class T>
void
quadgrid_t<T>::vtk_export (const char *filename,
			   const std::map<std::string, T> & f) const {

  std::ofstream ofs (filename, std::ofstream::out);

  // This is the XML format of a VTS file to write :

  ofs <<
    "<VTKFile type=\"StructuredGrid\" version=\"StructuredGrid\" byte_order=\"LittleEndian\">\n\
    <StructuredGrid WholeExtent=\"0 " << num_rows() << " 0 " << num_cols() << " 0 0\">\n \
      <Piece Extent=\"0 " << num_rows() << " 0 " << num_cols() << " 0 0\">\n";

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
  for (idx_t ii = 0; ii <= num_cols(); ++ii) {
    ofs << "          ";
    for (idx_t jj = 0; jj <= num_rows(); ++jj) {
      ofs << std::setprecision(16) << hx()*ii << " " << hy()*jj << " 0 ";
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


template <class T>
void
quadgrid_t<T>::octave_ascii_export
(const char *filename,
 std::map<std::string, T> const & vars) const {

  std::ofstream os (filename, std::ofstream::out);
  
  os << "# name: p" << std::endl
     << "# type: matrix" << std::endl
     << "# rows: 2" << std::endl
     << "# columns: " << num_owned_nodes () << std::endl;

  for (idx_t jj = 0; jj < num_cols () + 1; ++jj) {
    for (idx_t ii = 0; ii < num_rows () + 1; ++ii) {
      os  << std::setprecision(16) << jj*hx() << " ";
    }
  }
  os << std::endl;

  for (idx_t jj = 0; jj < num_cols () + 1; ++jj) {
    for (idx_t ii = 0; ii < num_rows () + 1; ++ii) {
      os  << std::setprecision(16) << ii*hy() << " ";
    }
  }
  os << std::endl;
  
  os << "# name: t" << std::endl
     << "# type: matrix" << std::endl
     << "# rows: 4" << std::endl
     << "# columns: " << num_local_cells () << std::endl;
  for (auto inode = 0;
       inode < quadgrid_t<std::vector<double>>::cell_t::nodes_per_cell;
       ++inode) {
    for (auto icell = begin_cell_sweep ();
	 icell != end_cell_sweep (); ++icell) {
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



