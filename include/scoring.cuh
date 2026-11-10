#ifndef SCORING_CUH
#define SCORING_CUH

#include <cstddef>

#include <cuda_runtime.h>

// Frameshift-aware local (Smith-Waterman-style) DNA-to-protein alignment.
// Rows index DNA bases (N), columns index protein residues (M); all matrices
// are row-major N x M. sc/ins/del hold scores, t_* hold traceback pointers.

// Sets column 0, rows 0-3, and the del-matrix boundary before scoring.
void scoring_init(char* input_seq, char* ref_seq, int* sc_mat, int* ins_mat, int* del_mat,
                  int* t_sc_mat, int* t_ins_mat, int* t_del_mat, size_t N, size_t M);

// CPU fill of rows 4..N-1.
void scoring_cpu(const char* input_seq, const char* ref_seq, int* sc_mat, int* ins_mat,
                 int* del_mat, int* t_sc_mat, int* t_ins_mat, int* t_del_mat, size_t N, size_t M);

// GPU fill of one anti-diagonal of tiles; launched once per currDiag.
__global__ void scoring_cuda(const char* __restrict__ input_seq, const char* __restrict__ ref_seq,
                             int* __restrict__ u_sc_mat, int* __restrict__ u_ins_mat,
                             int* __restrict__ u_del_mat, int* __restrict__ u_t_sc_mat,
                             int* __restrict__ u_t_ins_mat, int* __restrict__ u_t_del_mat,
                             unsigned int N, unsigned int M, unsigned int currDiag,
                             unsigned int submatrixSide, unsigned int currDiagFirstSubY);

#endif // SCORING_CUH
