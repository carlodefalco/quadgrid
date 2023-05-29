#include <json.hpp>
#include <particles.h>

#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <boost/iostreams/filter/gzip.hpp>
//#include <boost/iostreams/filter/bzip2.hpp>

int
main () {

  const std::string filename = "esempio.json";
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


  {
    std::map<std::string, std::vector<double>> vars
    {{"m", std::vector<double>(grid.num_global_nodes (), 0.)},
     {"vx", std::vector<double>(grid.num_global_nodes (), 0.)},
     {"vy", std::vector<double>(grid.num_global_nodes (), 0.)}
    };

    boost::iostreams::filtering_ostream outbuf;
    outbuf.push (boost::iostreams::gzip_compressor());
    // outbuf.push (boost::iostreams::bzip2_compressor());
    outbuf.push (boost::iostreams::file_sink (filename + fileext));
		
    nlohmann::json j(ptcls);
    j["grid_vars"] = vars;
    outbuf << std::setw(2) << j << std::endl;

    boost::iostreams::close (outbuf);
  }


  nlohmann::json j;
  boost::iostreams::filtering_istream inbuf;

  inbuf.push(boost::iostreams::gzip_decompressor());
  // inbuf.push(boost::iostreams::bzip2_decompressor());
  inbuf.push (boost::iostreams::file_source (filename + fileext));

  inbuf >> j;
  quadgrid_t<std::vector<double>> qg (j["grid_properties"]);
  particles_t p (j, qg);
  std::map<std::string, std::vector<double>> vars= 
    j["grid_vars"].get<std::map<std::string, std::vector<double>>> ();
    
  boost::iostreams::close (inbuf);
  
  {
    std::ofstream jf ("esempio2.json");
    nlohmann::json j (p);
    j["grid_vars"] = vars;
    jf << std::setw(2) << j;
    jf.close ();
  }
   
  return 0;
}
