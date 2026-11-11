#ifndef SCORING_CUH
#define SCORING_CUH

#include <cstddef>
#include <cstdint>

#include <cuda_runtime.h>

enum Trace : std::int8_t {
    TRACE_INS = -2,
    TRACE_DEL = -1,
    TRACE_STOP = 0,
    TRACE_FIRST = 1,
    TRACE_SHIFT_2 = 2,
    TRACE_MATCH = 3,
    TRACE_SHIFT_4 = 4,
};

void scoring_init(char* input_seq, char* ref_seq, int* sc_mat, int* ins_mat, int* del_mat,
                  int* t_sc_mat, int* t_ins_mat, int* t_del_mat, size_t N, size_t M);

void scoring_cpu(const char* input_seq, const char* ref_seq, int* sc_mat, int* ins_mat,
                 int* del_mat, int* t_sc_mat, int* t_ins_mat, int* t_del_mat, size_t N, size_t M);

__global__ void scoring_cuda(const char* __restrict__ codon_aa, const char* __restrict__ ref_seq,
                             int* __restrict__ u_sc_mat, int* __restrict__ u_ins_mat,
                             int* __restrict__ u_del_mat, int* __restrict__ u_t_sc_mat,
                             int* __restrict__ u_t_ins_mat, int* __restrict__ u_t_del_mat,
                             unsigned int N, unsigned int M, unsigned int currDiag,
                             unsigned int submatrixSide, unsigned int currDiagFirstSubY);

void upload_blosum62_tables();

#endif // SCORING_CUH
