#ifndef BLOSUM62_CUH
#define BLOSUM62_CUH

#include <cuda_runtime.h>

extern int blosum62mat[24][24];
extern __constant__ int d_blosum62mat[24][24];

extern int prot_to_idx[128];
extern __constant__ int d_prot_to_idx[128];

void initialize_protein_to_int_map();
void readBlosum62();

__host__ __device__ inline int place(const char a) {
#ifdef __CUDA_ARCH__
    const int* table = d_prot_to_idx;
#else
    const int* table = prot_to_idx;
#endif
    return (a >= 0 && a < 128) ? table[static_cast<unsigned char>(a)] : -1;
}

__host__ __device__ inline int score(const char a, const char b) {
#ifdef __CUDA_ARCH__
    return d_blosum62mat[place(a)][place(b)];
#else
    return blosum62mat[place(a)][place(b)];
#endif
}

#endif // BLOSUM62_CUH
