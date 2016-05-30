#include <cuda_runtime.h>

#include <iostream>

__global__ void init_phi_kernel(double *fab,
                                const int lo1, const int lo2, const int lo3,
                                const int hi1, const int hi2, const int hi3,
                                const double problo1, const double problo2, const double problo3,
                                const double probhi1, const double probhi2, const double probhi3,
                                const int jStride, const int kStride,
                                const int Nghost,
                                const double dx1, const double dx2, const double dx3);

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
  const double probhi1 = problo[0];
  const double probhi2 = problo[1];
  const double probhi3 = problo[2];

  const double dx1 = dx[0];
  const double dx2 = dx[1];
  const double dx3 = dx[2];

  // CUDA kernels can copy parameters by value to the GPU stack if and only if
  // they are primitive types. Anything more sophisticated (including arrays
  // and pointers to anything) must be managed using the runtime API.

  init_phi_kernel <<< 1, 1 >>> (fab,
                                lo1, lo2, lo3,
                                hi1, hi2, hi3,
                                problo1, problo2, problo3,
                                probhi1, probhi2, probhi3,
                                jStride, kStride,
                                Nghost,
                                dx1, dx2, dx3);

//  cudaDeviceSynchronize();

//  int i, j, k;
//  for (k = lo[2]; k <= hi[2]; ++k) {
//    for (j = lo[1]; j <= hi[1]; ++j) {
//      for (i = lo[0]; i <= hi[0]; ++i) {
//
//        const int ijk = (i+Nghost) + (j+Nghost)*jStride + (k+Nghost)*kStride;
//        fab[ijk] = 42.0;
//
//      }
//    }
//  }


  std::cout  << "done." << std::endl;
}
