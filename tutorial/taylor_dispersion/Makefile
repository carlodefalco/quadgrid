CPPFLAGS ?= -I../../include -I../json/single_include/nlohmann/ -I../json/single_include -I../cdf_timer/src
CXXFLAGS ?= -std=c++17 -O3
THRUST_CPU_CXXFLAGS ?= -std=c++20 -O3
QUADGRIDSRC ?= ../../src/particles.cpp
TIMERSRC ?= ../cdf_timer/src/timer.cpp
LDFLAGS ?=
LIBS ?= 


taylor_dispersion : taylor_dispersion.cpp taylor_dispersion.h
	nvcc $(CXXFLAGS) $(CPPFLAGS) -arch=sm_80 --expt-relaxed-constexpr -DUSE_THRUST -x cu $(QUADGRIDSRC) $(TIMERSRC)  taylor_dispersion.cpp -o taylordispersion $(LDFLAGS) $(LIBS) 

hip : taylor_dispersion.cpp taylor_dispersion.h
	hipcc $(CXXFLAGS) $(CPPFLAGS) -DUSE_THRUST $(QUADGRIDSRC) $(TIMERSRC) taylor_dispersion.cpp -o taylordispersion.hip $(LDFLAGS) $(LIBS)

thrust_cpu : taylor_dispersion.cpp taylor_dispersion.h
	nvcc $(THRUST_CPU_CXXFLAGS) $(CPPFLAGS) --expt-relaxed-constexpr -DUSE_THRUST -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP -DTHRUST_CPU $(QUADGRIDSRC) $(TIMERSRC) taylor_dispersion.cpp -o taylordispersion.cpu $(LDFLAGS) $(LIBS)

clean :
	$(RM) *.o *.csv *.vts *.json taylordispersion
