#include <cuda_runtime.h>

#include <CUDArena.H>

void*
CUDArena::alloc (std::size_t _sz)
{
    void *pt;
    cudaMallocManaged(&pt, _sz);
    cudaDeviceSynchronize();
    return pt;
}

void
CUDArena::free (void* pt)
{
    cudaDeviceSynchronize();
    cudaFree(pt);
}
