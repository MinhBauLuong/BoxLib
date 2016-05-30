__global__ void init_phi_kernel(double *fab, const int *lo, const int *hi, const int jStride, const int kStride, const int Nghost) {

  int i, j, k, ijk;

  for (k = lo[2]; k <= hi[2]; ++k) {
    for (j = lo[1]; j <= hi[1]; ++j) {
      for (i = lo[0]; i <= hi[0]; ++i) {

        ijk = (i+Nghost) + (j+Nghost)*jStride + (k+Nghost)*kStride;
        fab[ijk] = 42.0;

      }
    }
  }
}
