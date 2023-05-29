#include <particles.h>

#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <boost/iostreams/filter/gzip.hpp>
//#include <boost/iostreams/filter/zlib.hpp>
//#include <boost/iostreams/filter/bzip2.hpp>

int
main () {

  const std::string filename = "esempio.csv";
  const std::string fileext = ".gz";
  //const std::string fileext = ".bz";
  
  using idx_t = quadgrid_t<std::vector<double>>::idx_t;
  quadgrid_t<std::vector<double>> grid;
  grid.set_sizes (40, 60, 1./40., 1./60.);

  constexpr idx_t num_particles = 1500;
  particles_t ptcls (num_particles, {"label"}, {"m","vx","vy"}, grid);
  ptcls.dprops["m"].assign (num_particles, (1.) / num_particles );
  ptcls.dprops["vx"].assign (num_particles, 1.);
  ptcls.dprops["vy"].assign (num_particles, 1.);

  idx_t ilabel = 0;
  std::iota (ptcls.iprops["label"].begin (),
	     ptcls.iprops["label"].end (), ilabel);


  boost::iostreams::filtering_ostream outbuf;
  outbuf.push (boost::iostreams::gzip_compressor());
  //outbuf.push (boost::iostreams::zlib_compressor());
  // outbuf.push (boost::iostreams::bzip2_compressor());
  outbuf.push (boost::iostreams::file_sink (filename + fileext));
		
  ptcls.print<particles_t::output_format::csv> (outbuf);

  boost::iostreams::close (outbuf);
  return 0;
}

