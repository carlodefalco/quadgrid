# Taylor–Aris Dispersion Test Case

Taylor–Aris dispersion describes the spreading of a solute within a fluid due to the coupled effects of advection and diffusion.

---

## Prerequisites

Before compiling, ensure you load all required dependencies from the parent directory:

```bash
cd .. && make loaddeps
```

### Required Dependencies & Tools
* **`nvcc`** (NVIDIA CUDA Compiler — required for NVIDIA GPU target)
* **`hipcc`** (AMD ROCm / HIP Compiler — required for AMD GPU target)
* **`octave`** (Required for pre-processing input data generation)

---

## Compilation Guide

This test case supports three different execution backends. Choose the build target appropriate for your hardware environment:

| Hardware Target | Make Command | Notes & Configuration |
| :--- | :--- | :--- |
| **NVIDIA GPU** | `make taylordispersion` | *Update `-arch=sm_80` in the `Makefile` to match your GPU's compute capability before running.* |
| **CPU** | `make thrust_cpu` | Uses the Thrust host CPU backend. |
| **AMD GPU** | `make hip` | Uses the AMD HIP compiler environment. |

---

## Execution Instructions

### 1. Generate Input Data
Run the Octave script to setup initial solute distributions and simulation parameters:
```bash
octave write_td.m
```

### 2. Run the Simulation
Execute the compiled binary:
```bash
./taylordispersion
```
