## EXMPLE PSTL via oneDPC++ FLAGS :
## CPPFLAGS = "-DONEDPL_USE_TBB_BACKEND -I/opt/local/libexec/onetbb/include -I/opt/local/include/  "
## CXXFLAGS = "-O3 "
## LDFLAGS = "-L/opt/local/libexec/onetbb/lib/"

CPPFLAGS ?= -I../../include -I../json/single_include/nlohmann/ -I../json/single_include -I../cdf_timer/src
CXXFLAGS ?= -std=c++17 -O3
QUADGRIDSRC ?= ../../src/particles.cpp
TIMERSRC ?= ../cdf_timer/src/timer.cpp
LDFLAGS ?=  -L/opt/local/lib
LIBS ?= 


taylor_dispersion : taylor_dispersion.cpp taylor_dispersion.h
	nvcc $(CXXFLAGS) $(CPPFLAGS) -arch=sm_80 --expt-relaxed-constexpr -DUSE_THRUST -x cu $(QUADGRIDSRC) $(TIMERSRC)  taylor_dispersion.cpp -o taylordispersion $(LDFLAGS) $(LIBS) 

clean :
	$(RM) *.o *.csv *.vts *.json taylordispersion
