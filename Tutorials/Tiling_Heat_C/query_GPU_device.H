#ifndef QUERY_GPU_DEVICE_H
#define QUERY_GPU_DEVICE_H

#include <cuda_runtime.h>

void query_GPU_device(cudaDeviceProp &prop, const int device_idx);

#endif
