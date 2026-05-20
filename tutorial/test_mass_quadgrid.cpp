#include <json.hpp>
#include <particles.h>

#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#include "counter.h"
#include <timer.h>

// Explicit instantiation, dovrebbe generare tutti i metodi
template class quadgrid_t<std::vector<double>>;




int main(){

  cdf::timer::timer_t timer;

  constexpr auto filename = "mass.json";
  std::cout << std::setprecision(16);            // cifre significative


  nlohmann::json j;
  std::ifstream inbuf (filename);
  inbuf >> j;

  // create grid from fields in a json object
  quadgrid_t<std::vector<double>> qg (j["grid_properties"]);

  // create particles from properties in the json object
  // and the above created grid
  particles_t p (j, qg);
  


 std::map<std::string, std::vector<double>> vars= 
    j["grid_vars"].get<std::map<std::string, std::vector<double>>> ();
    
  inbuf.close ();



 p.init_particle_mesh ();


 std::vector<double> & mp=p.dprops["mp"];
 double mass_part=0;
  for (size_t ip=0; ip<p.num_particles;++ip){
    mass_part+=mp[ip];
  }
std::cout<<"mass_part_tot:"<<mass_part<<std::endl;



timer.tic("p2g");
p.p2g (vars, {"mp"}, {"Mn"});
timer.toc("p2g");
double mass_dof=0;
std::vector<double> & mdof=vars["Mn"];
for( size_t in=0; in< mdof.size();++in){

    mass_dof+=mdof[in];

}
std::cout<<"mass_dof_tot: "<<mass_dof<<std::endl;


std::cout<<"Test 2"<<std::endl;
double Mi=0.5;
// Test2 quantity conservation in g2p, only for constant fields
vars.at("Mn").assign(mdof.size(),Mi);

mass_dof=0;
for( size_t in=0; in< mdof.size();++in){

    mass_dof+=mdof[in];

}
std::cout<<"mass_dof_tot: "<<mass_dof<<std::endl;


mp.assign(p.num_particles,0.);
timer.tic("g2p");
p.g2p(vars, {"Mn"}, {"mp"});
timer.toc("g2p");

 mass_part=0;
  for (size_t ip=0; ip<p.num_particles;++ip){
    mass_part+=mp[ip];
  }

std::cout<<"mass_part_tot:"<<mass_part<<", expected: "<<Mi*p.num_particles<<std::endl;
std::cout<<"mdof size: "<<mdof.size()<<" num of dof: "<<qg.num_global_nodes()<<std::endl;
std::cout<<"Num particles: "<<p.num_particles<<std::endl;



mdof.assign(mdof.size(),0.);
timer.tic("p2g");
p.p2g(vars, {"mp"}, {"Mn"});
timer.toc("p2g");

mass_dof=0;
for( size_t in=0; in< mdof.size();++in){
    mass_dof+=mdof[in];

}
std::cout<<"mass_dof_tot: "<<mass_dof<<std::endl;



}















