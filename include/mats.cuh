#ifndef MATS_CUH
#define MATS_CUH

#include <cstddef>
#include <memory>

#include <cuda_runtime.h>

#include "cuda_check.cuh"

struct HostMats {
    size_t N, M;
    std::unique_ptr<int[]> sc, ins, del, t_sc, t_ins, t_del;

    HostMats(size_t N, size_t M)
        : N(N), M(M), sc(new int[N * M]), ins(new int[N * M]), del(new int[N * M]),
          t_sc(new int[N * M]), t_ins(new int[N * M]), t_del(new int[N * M]) {}

    size_t size() const {
        return N * M;
    }
};

struct DeviceMats {
    int *sc, *ins, *del, *t_sc, *t_ins, *t_del;
    unsigned int N, M;

    explicit DeviceMats(const HostMats& h)
        : N(static_cast<unsigned int>(h.N)), M(static_cast<unsigned int>(h.M)) {
        const size_t n = h.size();
        int** const dst[] = {&sc, &ins, &del, &t_sc, &t_ins, &t_del};
        const int* src[] = {h.sc.get(),   h.ins.get(),   h.del.get(),
                            h.t_sc.get(), h.t_ins.get(), h.t_del.get()};
        for (int k = 0; k < 6; k++) {
            checkCudaErrors(cudaMalloc(dst[k], n * sizeof(int)));
            checkCudaErrors(cudaMemcpy(*dst[k], src[k], n * sizeof(int), cudaMemcpyHostToDevice));
        }
    }
    void free() const {
        for (int* d : {sc, ins, del, t_sc, t_ins, t_del})
            checkCudaErrors(cudaFree(d));
    }
};

#endif // MATS_CUH
