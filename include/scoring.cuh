#ifndef SCORING_CUH
#define SCORING_CUH

#include <cstddef>
#include <cstdint>
#include <string>

#include <cuda_runtime.h>

#include "mats.cuh"

enum Trace : std::int8_t {
    TRACE_INS = -2,
    TRACE_DEL = -1,
    TRACE_STOP = 0,
    TRACE_FIRST = 1,
    TRACE_SHIFT_2 = 2,
    TRACE_MATCH = 3,
    TRACE_SHIFT_4 = 4,
};

void scoring_init(const std::string& input_seq, const std::string& ref_seq, HostMats& m);

void scoring_cpu(const std::string& input_seq, const std::string& ref_seq, HostMats& m);

__global__ void scoring_cuda(const char* __restrict__ codon_aa, const char* __restrict__ ref_seq,
                             DeviceMats m, unsigned int currDiag, unsigned int submatrixSide,
                             unsigned int currDiagFirstSubY);

void upload_blosum62_tables();

#endif // SCORING_CUH
