#ifndef CUDA_CHECK_CUH
#define CUDA_CHECK_CUH

#include <cstdio>
#include <cstdlib>

#include <cuda_runtime.h>

#define checkCudaErrors(val) check((val), #val, __FILE__, __LINE__)
template <typename T>
void check(T err, const char* const func, const char* const file, const int line) {
    if (err != cudaSuccess) {
        (void)fprintf(stderr, "CUDA error at: %s : %d\n", file, line);
        (void)fprintf(stderr, "%s %s\n", cudaGetErrorString(err), func);
        exit(1);
    }
}

#endif // CUDA_CHECK_CUH
