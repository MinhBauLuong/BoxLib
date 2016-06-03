#include <iostream>
#include <cuda_runtime.h>

void query_GPU_device() {
  int nDevices;

  cudaGetDeviceCount(&nDevices);

  std::cout << std::endl;
  std::cout << "Querying available CUDA devices ... " << std::endl << std::endl;
  for (int i = 0; i < nDevices; i++) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    std::cout << "Device Number: " << i << std::endl;
    std::cout << "  Device name: " << prop.name << std::endl;
    std::cout << "  Memory Clock Rate (MHz): " << int(double(prop.memoryClockRate)/1.0e3) << std::endl;
    std::cout << "  Memory Bus Width (bits): " << prop.memoryBusWidth << std::endl;
    std::cout << "  Peak Memory Bandwidth (GB/s): " << 2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6 << std::endl;
    std::cout << "  Size of L2 cache (KB) " << int(double(prop.l2CacheSize)/1024.0) << std::endl;
    std::cout << "  Supports managed memory? " << (prop.managedMemory == 0 ? "NO" : "YES") << std::endl;
    std::cout << "  Max size of each dimension of thread block: (" << prop.maxThreadsDim[0] << ", " << prop.maxThreadsDim[1] << ", " << prop.maxThreadsDim[2] << ")" << std::endl;
    std::cout << "  Max # of threads per block: " << prop.maxThreadsPerBlock << std::endl;
    std::cout << "  Max # of threads per multiprocessor: " << prop.maxThreadsPerMultiProcessor << std::endl;
    std::cout << "  Shared memory per block (KB): " << int(double(prop.sharedMemPerBlock)/1024.0) << std::endl;
    std::cout << "  Shared memory per multiprocessor (KB): " << int(double(prop.sharedMemPerMultiprocessor)/1024.0) << std::endl;
    std::cout << "  Global memory available on device (MB): " << int(double(prop.totalGlobalMem)/1024.0/1024.0) << std::endl;
    std::cout << "  Supports unified address space with host? " << (prop.unifiedAddressing == 0 ? "NO" : "YES") << std::endl;
    std::cout << "  Warp size (threads): " << prop.warpSize << std::endl;
    std::cout << std::endl;
  }
}
