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

#include <octave/oct.h>
#include <octave/interpreter.h>

#include <quadgrid_cpp.h>

#include <iostream>
#include <memory>

class
quadgrid
  : public octave_base_value {

public:

  quadgrid (octave_idx_type nx = 2, double hx = 1.0,
            octave_idx_type ny = 2, double hy = 1.0)
    : octave_base_value () {
    rep = std::make_shared<quadgrid_t<ColumnVector>> ();
    rep->set_sizes (ny, nx, hx, hy);
  };
  
  void
  print (std::ostream& os, bool pr_as_read_syntax = false) {
    os << "quadgrid object" << std::endl;
    os << "nx = " << rep->num_cols () << " ny = " << rep->num_rows () << std::endl;
    os << "hx = " << rep->hx () << " hy = " << rep->hy () << std::endl;
  }

  ~quadgrid () = default;

  bool
  is_defined (void) const
  { return true; }

  const std::shared_ptr<quadgrid_t<ColumnVector>>
  quadgrid_value (bool = false) const
  { return rep; }

  std::shared_ptr<quadgrid_t<ColumnVector>>
  quadgrid_value (bool = false) 
  { return rep; }
  
private:

  std::shared_ptr<quadgrid_t<ColumnVector>> rep;
  
  DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA;
};

const std::shared_ptr<quadgrid_t<ColumnVector>>
ov_quadgrid (const octave_value &in);

std::shared_ptr<quadgrid_t<ColumnVector>>
ov_quadgrid (octave_value &in);

bool
is_quadgrid (const octave_value &in);
