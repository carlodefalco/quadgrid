
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


dd : drift_diffusion.cpp drift_diffusion.h
	nvcc $(CXXFLAGS) $(CPPFLAGS) -DJSON_HAS_RANGES=0 -arch=sm_80 --expt-relaxed-constexpr -DUSE_THRUST -x cu $(QUADGRIDSRC) $(TIMERSRC)  drift_diffusion.cpp -o drift_diffusion $(LDFLAGS) $(LIBS) 

dd_cpu : drift_diffusion.cpp drift_diffusion.h
	nvcc $(CXXFLAGS)  $(CPPFLAGS) --expt-relaxed-constexpr -DTHRUST_CPU -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP -DUSE_THRUST  $(QUADGRIDSRC) $(TIMERSRC) drift_diffusion.cpp -o drift_diffusion $(LDFLAGS) $(LIBS)


clean :
	$(RM) *.o *.csv *.vts *.json drift_diffusion
