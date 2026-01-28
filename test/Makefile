# Makefile for quadgrid examples

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -I../include
LDFLAGS = 

# Directories
SRCDIR = ../src
INCLUDEDIR = ../include

# Source files from src/ that need to be compiled
SRC_OBJECTS = $(SRCDIR)/particles.o

# Target executable (change this for different examples)
MAIN = particle_sort_example

# Main source file (change this for different examples)
MAIN_SRC = particle_sort_example.cpp

# Object file for main source
MAIN_OBJ = $(MAIN_SRC:.cpp=.o)

# Default target
all: $(MAIN)

# Link the executable
$(MAIN): $(MAIN_OBJ) $(SRC_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compile main source file
$(MAIN_OBJ): $(MAIN_SRC)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Compile particles.cpp from src/
$(SRCDIR)/particles.o: $(SRCDIR)/particles.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean object files only
clean:
	rm *.o $(SRC_OBJECTS)

# Clean everything (object files and executables)
distclean: clean
	rm $(MAIN)

# Rebuild everything from scratch
rebuild: distclean all

.PHONY: all clean distclean rebuild
