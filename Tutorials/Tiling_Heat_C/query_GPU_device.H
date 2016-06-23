#ifndef QUERY_GPU_DEVICE_H
#define QUERY_GPU_DEVICE_H

#ifdef BL_CUDA

#include <cuda_runtime.h>

void query_GPU_device(cudaDeviceProp &prop, const int device_idx);

#endif /* BL_CUDA */

#endif
