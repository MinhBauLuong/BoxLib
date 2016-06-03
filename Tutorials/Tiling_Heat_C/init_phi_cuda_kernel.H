#ifndef INIT_PHI_CUDA_KERNEL_H
#define INIT_PHI_CUDA_KERNEL_H
__global__ void init_phi_kernel(double *fab,
                                const int lo1, const int lo2, const int lo3,
                                const int hi1, const int hi2, const int hi3,
                                const double problo1, const double problo2, const double problo3,
                                const double probhi1, const double probhi2, const double probhi3,
                                const int jStride, const int kStride,
                                const int Nghost,
                                const double dx1, const double dx2, const double dx3);
#endif
