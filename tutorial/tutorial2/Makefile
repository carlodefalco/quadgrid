
## OVERRIDE THE FLAGS THAT ARE NOT SUITABLE FOR YOUR MACHINE (DEFAULTS ARE DEFINED BELOW) :
## CPPFLAGS = 
## CXXFLAGS = 
## LDFLAGS  = 
## LIBS     =

## EXMPLE PSTL via oneDPC++ FLAGS :
## CPPFLAGS = "-DONEDPL_USE_TBB_BACKEND -I/opt/local/libexec/onetbb/include -I/opt/local/include/  "
## CXXFLAGS = "-O3 "
## LDFLAGS = "-L/opt/local/libexec/onetbb/lib/"


CXX ?= /opt/local/bin/clang++
CPPFLAGS ?= -I../../include -I../json/single_include/nlohmann/ -I../json/single_include -I../cdf_timer/src
CXXFLAGS ?= -std=c++20 -O3
QUADGRIDSRC ?= ../../src/particles.cpp
TIMERSRC ?= ../cdf_timer/src/timer.cpp
LDFLAGS ?=  -L/opt/local/lib
LIBS ?= 

thrust : grid_velocity_field_thrust.cpp
	nvcc -std=c++17 -O3 $(CPPFLAGS) -arch=sm_80 --expt-relaxed-constexpr -DUSE_THRUST -x cu $(QUADGRIDSRC) $(TIMERSRC)  grid_velocity_field_thrust.cpp -o grid_velocity_field $(LDFLAGS) $(LIBS) 

thrust_cpu : grid_velocity_field_thrust.cpp
	nvcc -std=c++20 -O3 $(CPPFLAGS) --expt-relaxed-constexpr -DTHRUST_CPU -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP -DUSE_THRUST  $(QUADGRIDSRC) $(TIMERSRC)  grid_velocity_field_thrust.cpp -o grid_velocity_field $(LDFLAGS) $(LIBS)

clean :
	$(RM) *.o *.csv *.vts *.json grid_velocity_field 
