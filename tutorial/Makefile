
## OVERRIDE THE FLAGS THAT ARE NOT SUITABLE FOR YOUR MACHINE (DEFAULTS ARE DEFINED BELOW) :
## CPPFLAGS = 
## CXXFLAGS = 
## LDFLAGS  = 
## LIBS     =

## EXMPLE PSTL via oneDPC++ FLAGS :
## CPPFLAGS = "-DONEDPL_USE_TBB_BACKEND -I/opt/local/libexec/onetbb/include -I/opt/local/include/  "
## CXXFLAGS = "-O3 "
## LDFLAGS = "-L/opt/local/libexec/onetbb/lib/"


CXX ?= /usr/bin/g++
CPPFLAGS ?= -I../include -I./json/single_include/nlohmann/ -I./json/single_include -I./cdf_timer/src
CXXFLAGS ?= -std=c++17 -O3
QUADGRIDSRC ?= ../src/particles.cpp
TIMERSRC ?= ./cdf_timer/src/timer.cpp
LDFLAGS ?=  -L/opt/local/lib
LIBS ?= 

grid_velocity_field : grid_velocity_field.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(QUADGRIDSRC) $(TIMERSRC) grid_velocity_field.cpp -o grid_velocity_field $(LDFLAGS) $(LIBS)

mass_test: test_mass_quadgrid.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(QUADGRIDSRC) $(TIMERSRC) test_mass_quadgrid.cpp -o mass_test
## clone dependencies via git
loaddeps :
	git clone --depth 1 git@github.com:carlodefalco/cdf_timer.git
	git clone --depth 1 git@github.com:nlohmann/json.git

clean :
	$(RM) *.o *.csv *.vts *.json grid_velocity_field
