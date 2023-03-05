#include <algorithm>
#include <random>
#include <quadgrid_cpp.h>
#include <particles.h>
#include <map>
#include <fstream>
#include <iostream>

using idx_t = quadgrid_t<std::vector<double>>::idx_t;

static constexpr double radius = .1;
static constexpr double xc = .5;
static constexpr double yc = .5;
static constexpr idx_t initial_num_particles = 1000000;

static bool
inside (double x, double y) {
  return ((x-xc)*(x-xc) + (y-yc)*(y-yc) < radius*radius);
};

int
main (int argc, char *argv[]) {

  quadgrid_t<std::vector<double>> grid;
  grid.set_sizes (16, 16, 1./16., 1./16.);

  particles_t ptcls (initial_num_particles, {"label"}, {"m", "vx", "vy"}, grid);
  ptcls.dprops["m"].assign (ptcls.num_particles, 1. / static_cast<double>(ptcls.num_particles));

  idx_t ilabel = 0;
  std::iota (ptcls.iprops["label"].begin (), ptcls.iprops["label"].end (), ilabel);


  std::cout << "before removale np = " << ptcls.num_particles << std::endl;

  std::ofstream o1 ("before_removal.csv", std::ios::out);
  ptcls.print<particles_t::output_format::csv> (o1);
  o1.close ();
  o1.open  ("before_removal.octtxt", std::ios::out);
  ptcls.print<particles_t::output_format::octave_ascii> (o1);
  o1.close ();

  ptcls.remove_in_region (inside);
  std::cout << "after removale np = " << ptcls.num_particles << std::endl;
  
  std::ofstream o2 ("after_removal.csv", std::ios::out);
  ptcls.print<particles_t::output_format::csv> (o2);
  o2.close ();
  o2.open ("after_removal.octtxt", std::ios::out);
  ptcls.print<particles_t::output_format::octave_ascii> (o2);
  o2.close ();
  
  return 0;
};


