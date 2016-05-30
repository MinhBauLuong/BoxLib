#include <cuda_runtime.h>

#include <iostream>

__global__ void init_phi_kernel(double *fab,
                                const int lo1, const int lo2, const int lo3,
                                const int hi1, const int hi2, const int hi3,
                                int jStride, const int kStride,
                                const int Nghost);

void init_phi(double *fab, const int *lo, const int *hi, const int jStride, const int kStride, const int Nghost)
{
  std::cout << "(" << lo[0] << ", " << lo[1] << ", " << lo[2] << ") -> (" << hi[0] << ", " << hi[1] << ", " << hi[2] << ") ... ";

  const int lo1 = lo[0];
  const int lo2 = lo[1];
  const int lo3 = lo[2];
  const int hi1 = hi[0];
  const int hi2 = hi[1];
  const int hi3 = hi[2];

  // CUDA kernels can copy parameters by value to the GPU stack if and only if
  // they are primitive types. Anything more sophisticated (including arrays
  // and pointers to anything) must be managed using the runtime API.

  init_phi_kernel <<< 1, 1 >>> (fab, lo1, lo2, lo3, hi1, hi2, hi3, jStride, kStride, Nghost);

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
