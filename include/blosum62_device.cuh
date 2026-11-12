#ifndef BLOSUM62_DEVICE_CUH
#define BLOSUM62_DEVICE_CUH

#include <cuda_runtime.h>

#include "blosum62.h"

__constant__ int d_blosum62mat[24][24];
__constant__ int d_prot_to_idx[128];

__host__ __device__ inline int place(const char a) {
#ifdef __CUDA_ARCH__
    const int* table = d_prot_to_idx;
#else
    const int* table = prot_to_idx;
#endif
    return (a >= 0 && a < 128) ? table[static_cast<unsigned char>(a)] : table['X'];
}

__host__ __device__ inline int score(const char a, const char b) {
#ifdef __CUDA_ARCH__
    return d_blosum62mat[place(a)][place(b)];
#else
    return blosum62mat[place(a)][place(b)];
#endif
}

#endif // BLOSUM62_DEVICE_CUH
