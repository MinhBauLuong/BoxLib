#include <cmath>

__global__ void init_phi_kernel(double *fab,
                                const int lo1, const int lo2, const int lo3,
                                const int hi1, const int hi2, const int hi3,
                                const double problo1, const double problo2, const double problo3,
                                const double probhi1, const double probhi2, const double probhi3,
                                const int jStride, const int kStride,
                                const int Nghost,
                                const double dx1, const double dx2, const double dx3) {

  int i, j, k, ijk_fab;

  double x, y, z, r2;

  // Convert CUDA thread indices into indicies of the FAB that each thread will
  // modify. Isn't this beautiful??

  i = lo1 + (blockIdx.x * blockDim.x) + threadIdx.x;
  j = lo2 + (blockIdx.y * blockDim.y) + threadIdx.y;
  k = lo3 + (blockIdx.z * blockDim.z) + threadIdx.z;

  // WARNING: for now this assumes that the thread block and the FAB have
  // identical dimensions!!

  x = problo1 + (double(i)+0.5) * dx1;
  y = problo2 + (double(j)+0.5) * dx2;
  z = problo3 + (double(k)+0.5) * dx3;

  r2 = ((x-0.25)*(x-0.25) + (y-0.25)*(y-0.25) + (z-0.25)*(z-0.25)) * 100.0;

//  ijk_fab = (i+Nghost) + (j+Nghost)*jStride + (k+Nghost)*kStride;
  ijk_fab = i + j*jStride + k*kStride;

  fab[ijk_fab] = 1.0 + std::exp(-r2);

}
