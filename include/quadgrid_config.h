#ifndef QUADGRID_CONFIG_H
#define QUADGRID_CONFIG_H
#endif
#ifdef USE_MPI_H
#include <mpi.h>
#else
#define MPI_Comm int
#define MPI_COMM_WORLD 1
#define MPI_Initialized(x)
#define MPI_Comm_size(x, y) { *y = 0; }
#define MPI_Comm_rank(x, y) { *y = 0; }
#endif
#ifdef USE_THRUST
#include<thrust/host_vector.h>
#define DEVICE  __device__
#define HOST  __host__
#define vector_t  thrust::host_vector
#define algorithm_namespace thrust
#else
#include<atomicAdd.h>
#include<vector>
#define DEVICE  
#define HOST  
#define vector_t std::vector
#define algorithm_namespace std
#endif
using real_t = double;
