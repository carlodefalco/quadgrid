#include <mmorton.h>

static constexpr morton_code_t B[] = {0x5555555555555555, 0x3333333333333333, 0x0F0F0F0F0F0F0F0F, 0x00FF00FF00FF00FF, 0x0000FFFF0000FFFF};
static constexpr morton_code_t M[] = {0x4444444444444444, 0x3030303030303030, 0x0F000F000F000F00, 0x00FF000000FF0000, 0x0000FFFF00000000};
static constexpr morton_code_t S[] = {1, 2, 4, 8, 16};

morton_code_t
coord_2_morton (coord_t _x, coord_t _y) {

  morton_code_t x = static_cast<morton_code_t> (_x);
  morton_code_t y = static_cast<morton_code_t> (_y);
  
  x = (x | (x << S[3])) & B[3];
  x = (x | (x << S[2])) & B[2];
  x = (x | (x << S[1])) & B[1];
  x = (x | (x << S[0])) & B[0];

  y = (y | (y << S[3])) & B[3];
  y = (y | (y << S[2])) & B[2];
  y = (y | (y << S[1])) & B[1];
  y = (y | (y << S[0])) & B[0];

  return x | (y << 1);
}

// Below is a slower but easier to understand implementation
// the two versions of this function should return the same results
// which can be used for testing and debugging
/*
morton_code_t
coord_2_morton (coord_t x, coord_t y) {

  static coord_t W = 0;
  static coord_t H = 0;
  static morton_code_t xmask = 0;
  static morton_code_t ymask = 0;
  
  morton_code_t res = 0;
  
  W = (1 << (max_level -1));
  H = (1 << (max_level -1));

  xmask = 1 << (2 * (max_level - 1));
  ymask = xmask << 1;

while (x || y) {

    
    std::cout << "res = " << res << " x = " << x << " y = " << y << std::endl;
    std::cout << "H = " << H << " W = " << W << std::endl;
    std::cout << "xmask = " << xmask << " ymask = " << ymask << std::endl;
    
    
    if (y / H > 0) {
      res |= ymask;
      y ^= H;
    }
    
    if (x / W > 0) {
      res |= xmask;
      x ^= W;
    }

    xmask = xmask >> 2;
    ymask = ymask >> 2;
    H = H >> 1;
    W = W >> 1;
  }
  
  return (res);
}
*/

void
morton_2_coord (morton_code_t code, coord_t &_x, coord_t &_y) {

  morton_code_t x = code        & B[0];
  morton_code_t y = (code >> 1) & B[0];

  x = x | ((x & M[0]) >> S[0]);
  x = x | ((x & M[1]) >> S[1]);
  x = x | ((x & M[2]) >> S[2]);
  x = x | ((x & M[3]) >> S[3]);
  x = x | ((x & M[4]) >> S[4]);
  _x = x;
  
  y = y | ((y & M[0]) >> S[0]);
  y = y | ((y & M[1]) >> S[1]);
  y = y | ((y & M[2]) >> S[2]);
  y = y | ((y & M[3]) >> S[3]);
  y = y | ((y & M[4]) >> S[4]);
  _y = y;
}

// Below is a slower but easier to understand implementation
// the two versions of this function should return the same results
// which can be used for testing and debugging
/*
void
morton_2_coord (morton_code_t code, coord_t &x, coord_t &y) {

  static coord_t level  = 0;
  
  x = 0;
  y = 0;
  level = 0;
  
  while (code) {
    x |= (code & 1) << level;
    code = code >> 1;
    y |= (code & 1) << level++;
    code = code >> 1;
  }
  
}
*/

