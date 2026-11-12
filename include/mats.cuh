#ifndef MATS_CUH
#define MATS_CUH

#include <cstddef>
#include <memory>

#include <cuda_runtime.h>

#include "cuda_check.cuh"
#include "scoring.cuh"

struct HostMats {
    std::unique_ptr<int[]> sc, ins, del, t_sc, t_ins, t_del;

    explicit HostMats(size_t n)
        : sc(new int[n]), ins(new int[n]), del(new int[n]), t_sc(new int[n]), t_ins(new int[n]),
          t_del(new int[n]) {}

    void init(const char* dna, const char* prot, size_t N, size_t M) {
        scoring_init(dna, prot, sc.get(), ins.get(), del.get(), t_sc.get(), t_ins.get(),
                     t_del.get(), N, M);
    }
};

struct DeviceMats {
    int *sc, *ins, *del, *t_sc, *t_ins, *t_del;

    DeviceMats(const HostMats& h, size_t n) {
        int** const dst[] = {&sc, &ins, &del, &t_sc, &t_ins, &t_del};
        const int* src[] = {h.sc.get(),   h.ins.get(),   h.del.get(),
                            h.t_sc.get(), h.t_ins.get(), h.t_del.get()};
        for (int k = 0; k < 6; k++) {
            checkCudaErrors(cudaMalloc(dst[k], n * sizeof(int)));
            checkCudaErrors(cudaMemcpy(*dst[k], src[k], n * sizeof(int), cudaMemcpyHostToDevice));
        }
    }
    ~DeviceMats() {
        for (int* d : {sc, ins, del, t_sc, t_ins, t_del})
            checkCudaErrors(cudaFree(d));
    }
    DeviceMats(const DeviceMats&) = delete;
    DeviceMats& operator=(const DeviceMats&) = delete;
};

#endif // MATS_CUH
