quadgrid - a simple c++ library for particles in a cartesian quad grid
======

This library provides simple (template) classes that for a useful
starting point for the implementation of [MPM/PIC](https://www.sciencedirect.com/science/article/abs/pii/S0065215620300120) methods 


The repository consists of the following folders :

* `include`  contains headers (*`.h`) and template method definition
  files (`*_imp.h`) the main files included are
	* `quadgrid_cpp.h`  including the declaration of the (template) class
      `quadgrid_t` representing the quad grid
	* `quadgrid_cpp_imp.h` contains *out-of-line* definitions for
      template methods of the `quadgrid_t`  class   
	* `particles.h` declares the `particles_t` clares representing
      particles embedded in a `quadgrid_t` grid
* `src` contains implementation of methods in the above classes that
  do not depend on template parameters
* `test`  provides a few tests and examples
* `octave` provides a draft of an interface for accessing qudtree
  objects from within the [GNU Octave](http://www.octave.org)
  interpreter, which consists of 
    * `quadgrid.h`  defining the `quadgrid` class inheriting from
      `octave_base_value`
	* `quadgrid.cc` defines two Octave functions `quadgrid` and
      `quadgrid_loop`accessible from the interpreter

### Building the examples

To build the examples move to the `test` directory and run

    mpicxx -std=c++17 -I../include -o particle_sort_example particle_sort_example.cpp ../src/particles.cpp
    
### Main methods in the particles_t class

* `particles_t::p2g` implements transfer of quantities from the
  particles to the grid according to the formula
  
$$
  u_i = \sum_p N_i(x_p) U_p 
 $$


* `particles_t::g2p` implements transfer of quantities from the
  particles to the grid according to the formula
  
$$
 U_p = \sum_i N_i(x_p) u_i  
 $$
 
 
