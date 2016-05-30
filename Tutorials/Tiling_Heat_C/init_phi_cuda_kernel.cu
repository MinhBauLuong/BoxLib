__global__ void init_phi_kernel(double *fab,
                                const int lo1, const int lo2, const int lo3,
                                const int hi1, const int hi2, const int hi3,
                                const int jStride, const int kStride,
                                const int Nghost) {

  int i, j, k, ijk;

  for (k = lo3; k <= hi3; ++k) {
    for (j = lo2; j <= hi2; ++j) {
      for (i = lo1; i <= hi1; ++i) {

        ijk = (i+Nghost) + (j+Nghost)*jStride + (k+Nghost)*kStride;
        fab[ijk] = 42.0;

      }
    }
  }

}
