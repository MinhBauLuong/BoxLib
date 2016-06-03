#include <cmath>

__global__ void init_phi_kernel(double *fab,
                                const int lo1, const int lo2, const int lo3,
                                const int hi1, const int hi2, const int hi3,
                                const double problo1, const double problo2, const double problo3,
                                const double probhi1, const double probhi2, const double probhi3,
                                const int jStride, const int kStride,
                                const int Nghost,
                                const double dx1, const double dx2, const double dx3) {

  int i, j, k, ijk;

  double x, y, z, r2;

  for (k = lo3; k <= hi3; ++k) {
    z = problo3 + (double(k)+0.5) * dx3;
    for (j = lo2; j <= hi2; ++j) {
      y = problo2 + (double(j)+0.5) * dx2;
      for (i = lo1; i <= hi1; ++i) {
        x = problo1 + (double(i)+0.5) * dx1;

        r2 = ((x-0.25)*(x-0.25) + (y-0.25)*(y-0.25) + (z-0.25)*(z-0.25)) * 100.0;

        ijk = (i+Nghost) + (j+Nghost)*jStride + (k+Nghost)*kStride;

        fab[ijk] = 1.0 + std::exp(-r2);

      }
    }
  }

}
