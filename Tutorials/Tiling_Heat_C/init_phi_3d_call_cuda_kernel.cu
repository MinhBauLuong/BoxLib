#include <cuda_runtime.h>

#include <iostream>

__global__ void init_phi_kernel(double *fab, const int *lo, const int *hi, const int jStride, const int kStride, const int Nghost);

void init_phi(double *fab, const int *lo, const int *hi, const int jStride, const int kStride, const int Nghost)
{
  std::cout << "(" << lo[0] << ", " << lo[1] << ", " << lo[2] << ") -> (" << hi[0] << ", " << hi[1] << ", " << hi[2] << ") ... ";

  init_phi_kernel<<<1,1>>>(fab,lo,hi,jStride,kStride,Nghost);

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
