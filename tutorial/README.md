# QuadGrid Tutorials

This directory contains example problems and tutorials demonstrating
the usage of 
**QuadGrid**—a C++/Octave library for simple Cartesian quad grids. 
These tutorials cover advection-diffusion problems on both CPU and GPU (via Thrust).


## Overview

### Tutorial 1
Solves a 2D transport problem combining deterministic velocity advection with stochastic Brownian motion. Backend: CPU.

### Tutorial 2
This tutorial demonstrates the use of the **Thrust** library to solve a 2D advection only problem. Backend CPU(Standard/Thrust), NVIDIA GPU(Thrust).

### Tutorial 3 
Solves a 2D transport problem combining deterministic velocity advection and deterministic diffusion. Backend: CPU(Thrust), NVIDIA GPU(Thrust). 

### Tutorial 4: Taylor Dispersion
Simulates solute spreading in shear flows. Reference: *https://doi.org/10.1063/1.3078518*. Backend: CPU(Thrust), NVIDIA GPU(Thrust), AMD GPU(Thrust + hip).

