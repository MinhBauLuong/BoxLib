#ifndef INIT_PHI_3D_CALL_CUDA_KERNEL_H
#define INIT_PHI_3D_CALL_CUDA_KERNEL_H
void init_phi(double *fab,
              const int *lo, const int *hi,
              const double *problo, const double *probhi,
              const int jStride, const int kStride,
              const int Nghost,
              const double *dx);
#endif
