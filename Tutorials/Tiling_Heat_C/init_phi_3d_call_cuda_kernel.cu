#include <iostream>

#include <cuda_runtime.h>

#include <init_phi_cuda_kernel.H>

void init_phi(double *fab,
              const int *lo, const int *hi,
              const double *problo, const double *probhi,
              const int jStride, const int kStride,
              const int Nghost,
              const double *dx)
{
  std::cout << "(" << lo[0] << ", " << lo[1] << ", " << lo[2] << ") -> (" << hi[0] << ", " << hi[1] << ", " << hi[2] << ") ... ";

  const int lo1 = lo[0];
  const int lo2 = lo[1];
  const int lo3 = lo[2];
  const int hi1 = hi[0];
  const int hi2 = hi[1];
  const int hi3 = hi[2];

  const double problo1 = problo[0];
  const double problo2 = problo[1];
  const double problo3 = problo[2];
  const double probhi1 = probhi[0];
  const double probhi2 = probhi[1];
  const double probhi3 = probhi[2];

  const double dx1 = dx[0];
  const double dx2 = dx[1];
  const double dx3 = dx[2];

  // CUDA kernels can copy parameters by value to the GPU stack if and only if
  // they are primitive types. Anything more sophisticated (including arrays
  // and pointers to anything) must be managed using the runtime API.

  // This configuration guarantees that we can operate on exactly a 64^3 box.
  dim3 threadsPerBlock(8,8,8);
  dim3 blocksInGrid(8,8,8);

  init_phi_kernel <<< blocksInGrid, threadsPerBlock >>> (fab,
                                lo1, lo2, lo3,
                                hi1, hi2, hi3,
                                problo1, problo2, problo3,
                                probhi1, probhi2, probhi3,
                                jStride, kStride,
                                Nghost,
                                dx1, dx2, dx3);

  std::cout  << "done." << std::endl;
}
