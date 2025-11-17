#ifndef HAVE_MMORTON_H
#define HAVE_MMORTON_H

#include <cstdint>
#include <numeric>

using coord_t               = uint_least32_t;
using morton_code_t         = uint_fast64_t;

constexpr coord_t max_level = 31;
constexpr coord_t max_coord = std::numeric_limits<coord_t>::max ();

morton_code_t
coord_2_morton (coord_t _x, coord_t _y); 

void
morton_2_coord (morton_code_t code, coord_t &_x, coord_t &_y); 

#endif // HAVE_MMORTON_H
