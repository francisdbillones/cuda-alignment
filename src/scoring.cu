#include "scoring.cuh"

#ifdef __INTELLISENSE__
#define __CUDACC__
#endif // __INTELLISENSE__
#include <cooperative_groups.h>
#ifdef __INTELLISENSE__
#undef __CUDACC__
#endif // __INTELLISENSE__

#include "blosum62_device.cuh"
#include "codon.cuh"
#include "cuda_check.cuh"

namespace cg = cooperative_groups;

constexpr int gep = 2;   // extend penalty
constexpr int gop = 3;   // opening penalty
constexpr int shift = 4; // shift penalty
constexpr int infn = -999;

__constant__ int c_gep = 2;   // extend penalty
__constant__ int c_gop = 3;   // opening penalty
__constant__ int c_shift = 4; // shift penalty

void upload_blosum62_tables() {
    checkCudaErrors(cudaMemcpyToSymbol(d_blosum62mat, blosum62mat, sizeof(blosum62mat), 0,
                                       cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpyToSymbol(d_prot_to_idx, prot_to_idx, sizeof(prot_to_idx), 0,
                                       cudaMemcpyHostToDevice));
}

void scoring_init(char* input_seq, char* ref_seq, int* sc_mat, int* ins_mat, int* del_mat,
                  int* t_sc_mat, int* t_ins_mat, int* t_del_mat, size_t N, size_t M) {
    int insert = 0;
    int del = 0;
    int xscore = 0;
    char prot_seq;

    for (size_t i = 0; i < N; i++) {
        ins_mat[i * M + 0] = infn;
        t_ins_mat[i * M + 0] = infn;
    }

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            if (i == 0 || j == 0) {
                sc_mat[i * M + j] = 0;
                t_sc_mat[i * M + j] = 0;
            }
        }
    }

    for (size_t j = 0; j < M; j++) {
        del_mat[0 * M + j] = infn;
        del_mat[2 * M + j] = infn;
        del_mat[3 * M + j] = infn;
        del_mat[1 * M + j] = sc_mat[0 * M + j] - gop - gep;

        t_del_mat[0 * M + j] = infn;
        t_del_mat[2 * M + j] = infn;
        t_del_mat[3 * M + j] = infn;
        t_del_mat[1 * M + j] = 1;
    }

    for (size_t i = 4; i < N; i++) {
        del = del_mat[(i - 3) * M + 0] - gep;
        xscore = sc_mat[(i - 3) * M + 0] - gop - gep;

        del_mat[i * M + 0] = (del > xscore) ? del : xscore;
    }

    for (size_t i = 0; i < 4; i++) {
        for (size_t j = 1; j < M; j++) {
            insert = ins_mat[i * M + (j - 1)] - gep;
            xscore = sc_mat[i * M + (j - 1)] - gop - gep;
            prot_seq = translate_codon(input_seq + i - 1);

            if (insert > xscore) {
                ins_mat[i * M + j] = insert;
                t_ins_mat[i * M + j] = 0;
            } else {
                ins_mat[i * M + j] = xscore;
                t_ins_mat[i * M + j] = 1;
            }

            insert = ins_mat[i * M + j];
            del = del_mat[i * M + j];

            if (i == 1) {
                xscore = sc_mat[0 * M + (j - 1)] + score(prot_seq, ref_seq[j - 1]);
                if (insert >= del && insert >= xscore) {
                    sc_mat[i * M + j] = insert;
                    t_sc_mat[i * M + j] = -2;
                } else if (del >= insert && del >= xscore) {
                    sc_mat[i * M + j] = del;
                    t_sc_mat[i * M + j] = -1;
                } else {
                    sc_mat[i * M + j] = xscore;
                    t_sc_mat[i * M + j] = 1;
                }
            } else if (i == 2) {
                xscore = sc_mat[0 * M + (j - 1)] + score(prot_seq, ref_seq[j - 1]) - shift;
                if (insert > xscore) {
                    sc_mat[i * M + j] = insert;
                    t_sc_mat[i * M + j] = -2;
                } else {
                    sc_mat[i * M + j] = xscore;
                    t_sc_mat[i * M + j] = 2;
                }
            } else if (i == 3) {
                xscore = sc_mat[1 * M + (j - 1)] + score(prot_seq, ref_seq[j - 1]) - shift;
                if (insert > xscore) {
                    sc_mat[i * M + j] = insert;
                    t_sc_mat[i * M + j] = -2;
                } else {
                    sc_mat[i * M + j] = xscore;
                    t_sc_mat[i * M + j] = 2;
                }
            }

            if (sc_mat[i * M + j] < 0)
                sc_mat[i * M + j] = 0;
        }
    }
}

void scoring_cpu(const char* input_seq, const char* ref_seq, int* sc_mat, int* ins_mat,
                 int* del_mat, int* t_sc_mat, int* t_ins_mat, int* t_del_mat, size_t N, size_t M) {
    int insert = 0;
    int del = 0;
    int xscore = 0;
    int sc_1 = 0, sc_2 = 0, sc_3 = 0;
    int scoring = 0;
    char prot_seq;

    for (int i = 4; i < N; i++) {
        prot_seq = translate_codon(input_seq + i - 1);
        for (int j = 1; j < M; j++) {
            scoring = score(prot_seq, ref_seq[j - 1]);
            insert = ins_mat[i * M + (j - 1)] - gep;
            xscore = sc_mat[i * M + (j - 1)] - gop - gep;

            ins_mat[i * M + j] = (insert > xscore) ? insert : xscore;
            t_ins_mat[i * M + j] = (insert > xscore) ? 0 : 1;

            del = del_mat[(i - 3) * M + j] - gep;
            xscore = sc_mat[(i - 3) * M + j] - gop - gep;

            del_mat[i * M + j] = (del > xscore) ? del : xscore;
            t_del_mat[i * M + j] = (del > xscore) ? 0 : 1;

            if (i < N - 1) {
                insert = ins_mat[i * M + j];
                del = del_mat[i * M + j];
                sc_1 = sc_mat[(i - 2) * M + (j - 1)] + scoring - shift;
                sc_2 = sc_mat[(i - 3) * M + (j - 1)] + scoring;
                sc_3 = sc_mat[(i - 4) * M + (j - 1)] + scoring - shift;

                if (insert >= del && insert >= sc_1 && insert >= sc_2 && insert >= sc_3) {
                    sc_mat[i * M + j] = insert;
                    t_sc_mat[i * M + j] = -2;
                } else if (del >= insert && del >= sc_1 && del >= sc_2 && del >= sc_3) {
                    sc_mat[i * M + j] = del;
                    t_sc_mat[i * M + j] = -1;
                } else if (sc_2 >= insert && sc_2 >= del && sc_2 >= sc_1 && sc_2 >= sc_3) {
                    sc_mat[i * M + j] = sc_2;
                    t_sc_mat[i * M + j] = 3;
                } else if (sc_1 >= insert && sc_1 >= del && sc_1 >= sc_2 && sc_1 >= sc_3) {
                    sc_mat[i * M + j] = sc_1;
                    t_sc_mat[i * M + j] = 2;
                } else if (sc_3 >= insert && sc_3 >= del && sc_3 >= sc_1 && sc_3 >= sc_2) {
                    sc_mat[i * M + j] = sc_3;
                    t_sc_mat[i * M + j] = 4;
                }
            }

            if (sc_mat[i * M + j] < 0) {
                sc_mat[i * M + j] = 0;
            }

            if (i == N - 1) {
                ins_mat[i * M + j] = infn;
                sc_mat[i * M + j] = 0;
                t_sc_mat[i * M + j] = infn;
            }
        }
    }
}

__global__ void scoring_cuda(const char* __restrict__ codon_aa, const char* __restrict__ ref_seq,
                             int* __restrict__ u_sc_mat, int* __restrict__ u_ins_mat,
                             int* __restrict__ u_del_mat, int* __restrict__ u_t_sc_mat,
                             int* __restrict__ u_t_ins_mat, int* __restrict__ u_t_del_mat,
                             unsigned int N, unsigned int M, unsigned int currDiag,
                             unsigned int submatrixSide, unsigned int currDiagFirstSubY) {
    cg::thread_block block = cg::this_thread_block();

    unsigned int diagBlockIdx = blockIdx.x;

    unsigned int submatrixStartY = currDiagFirstSubY + diagBlockIdx;
    unsigned int submatrixStartX = currDiag - submatrixStartY;

    int i = submatrixStartY * blockDim.y + threadIdx.y;
    int j = submatrixStartX * blockDim.x + threadIdx.x;

    int insert = 0;
    int del = 0;
    int xscore = 0;
    int sc_1 = 0, sc_2 = 0, sc_3 = 0;
    int scoring = 0;
    char prot_seq;

    if (i >= 4 && i < N && j >= 1 && j < M) {
        for (unsigned int sub_diag_in_tile = 0; sub_diag_in_tile < 2 * submatrixSide - 1;
             ++sub_diag_in_tile) {
            if (((i - (submatrixStartY * blockDim.y)) + (j - (submatrixStartX * blockDim.x))) ==
                sub_diag_in_tile) {
                if (i == N - 1) {
                    u_ins_mat[i * M + j] = infn;
                    u_sc_mat[i * M + j] = 0;
                    u_t_sc_mat[i * M + j] = infn;
                } else {
                    prot_seq = codon_aa[i];
                    scoring = score(prot_seq, ref_seq[j - 1]);
                    insert = u_ins_mat[i * M + (j - 1)] - c_gep;
                    xscore = u_sc_mat[i * M + (j - 1)] - c_gop - c_gep;

                    u_ins_mat[i * M + j] = (insert > xscore) ? insert : xscore;
                    u_t_ins_mat[i * M + j] = (insert > xscore) ? 0 : 1;

                    del = u_del_mat[(i - 3) * M + j] - c_gep;
                    xscore = u_sc_mat[(i - 3) * M + j] - c_gop - c_gep;

                    u_del_mat[i * M + j] = (del > xscore) ? del : xscore;
                    u_t_del_mat[i * M + j] = (del > xscore) ? 0 : 1;

                    insert = u_ins_mat[i * M + j];
                    del = u_del_mat[i * M + j];
                    sc_1 = u_sc_mat[(i - 2) * M + (j - 1)] + scoring - c_shift;
                    sc_2 = u_sc_mat[(i - 3) * M + (j - 1)] + scoring;
                    sc_3 = u_sc_mat[(i - 4) * M + (j - 1)] + scoring - c_shift;

                    int max_current_score;
                    int trace_value;

                    max_current_score = (insert >= del) ? insert : del;
                    trace_value = (insert >= del) ? -2 : -1;

                    if (sc_2 > max_current_score) {
                        max_current_score = sc_2;
                        trace_value = 3;
                    }
                    if (sc_1 > max_current_score) {
                        max_current_score = sc_1;
                        trace_value = 2;
                    }
                    if (sc_3 > max_current_score) {
                        max_current_score = sc_3;
                        trace_value = 4;
                    }

                    u_sc_mat[i * M + j] = (max_current_score < 0) ? 0 : max_current_score;
                    u_t_sc_mat[i * M + j] = trace_value;
                }
            }
            cg::sync(block);
        }
    }
}
