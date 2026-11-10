#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <algorithm>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#ifdef __INTELLISENSE__
#define __CUDACC__
#endif // __INTELLISENSE__
#include <cooperative_groups.h>
#ifdef __INTELLISENSE__
#undef __CUDACC__
#endif // __INTELLISENSE__

#include "blosum62.cuh"
#include "codon.cuh"
#include "cuda_check.cuh"
#include "fasta.h"
#include "gputimer.h"

using std::cerr;
using std::cin;
using std::copy;
using std::cout;
using std::endl;
using std::milli;
using std::min;
using std::string;
using std::to_string;
using std::vector;
using std::chrono::duration;
using std::chrono::steady_clock;

namespace cg = cooperative_groups;

constexpr int gep = 2;   // extend penalty
constexpr int gop = 3;   // opening penalty
constexpr int shift = 4; // shift penalty
constexpr int infn = -999;
string myArray[40000][5];

__constant__ int c_gep = 2;   // extend penalty
__constant__ int c_gop = 3;   // opening penalty
__constant__ int c_shift = 4; // shift penalty
int score_top[4];

int insert = 0;
int del = 0;
int xscore = 0;
int end = 3;
int sc_1 = 0, sc_2 = 0, sc_3 = 0;
int scoring = 0;
char prot_seq;

void init_local_v2(char* input_seq, char* ref_seq, int* sc_mat, int* ins_mat, int* del_mat,
                   int* t_sc_mat, int* t_ins_mat, int* t_del_mat, size_t N, size_t M) {
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

void scoring_local_v2(const char* input_seq, const char* ref_seq, int* sc_mat, int* ins_mat,
                      int* del_mat, int* t_sc_mat, int* t_ins_mat, int* t_del_mat, size_t N,
                      size_t M) {
    for (int i = 4; i < N; i++) {
        for (int j = 1; j < M; j++) {
            prot_seq = translate_codon(input_seq + i - 1);
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

__global__ void scoring_local_v2_cuda(const char* __restrict__ input_seq,
                                      const char* __restrict__ ref_seq, int* __restrict__ u_sc_mat,
                                      int* __restrict__ u_ins_mat, int* __restrict__ u_del_mat,
                                      int* __restrict__ u_t_sc_mat, int* __restrict__ u_t_ins_mat,
                                      int* __restrict__ u_t_del_mat, unsigned int N, unsigned int M,
                                      unsigned int currDiag, unsigned int submatrixSide,
                                      unsigned int currDiagFirstSubY) {
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
                    prot_seq = translate_codon(input_seq + i - 1);
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

void top5(int score, int index, int top_i, int top_j, int* score_top, int* top_i_max,
          int* top_j_max, int* top_indices) {
    for (int x = 0; x < 5; x++) {
        if (score >= score_top[x]) {
            for (int y = 4; y > x; y--) {
                score_top[y] = score_top[y - 1];
                top_i_max[y] = top_i_max[y - 1];
                top_j_max[y] = top_j_max[y - 1];
                top_indices[y] = top_indices[y - 1];
            }

            score_top[x] = score;
            top_i_max[x] = top_i;
            top_j_max[x] = top_j;
            top_indices[x] = index;

            break;
        }
    }
}

void routine(int trace, int& i, int& j, char* str, char* ref_seq, string& final_seq1,
             string& final_seq2, string& frameshift) {
    int k = 0;

    if (trace == 3) {
        i -= 3;
        j--;
        k = i + 3;
        final_seq1 += translate_codon(str + k - 1);
        final_seq2 += ref_seq[j];
        frameshift += " ";
    } else if (trace == -2) {
        j--;
        final_seq1 += "-";
        final_seq2 += ref_seq[j];
        frameshift += " ";
    } else if (trace == -1) {
        i -= 3;
        k = i + 3;
        final_seq1 += translate_codon(str + k - 1);
        final_seq2 += "-";
        frameshift += " ";
    } else if (trace == 4) {
        i -= 4;
        j--;
        k = i + 4;
        final_seq1 += translate_codon(str + k - 1);
        final_seq2 += ref_seq[j];
        frameshift += "*";

    } else if (trace == 2) {
        i -= 2;
        j--;
        k = i + 2;
        final_seq1 += translate_codon(str + k - 1);
        final_seq2 += ref_seq[j];
        frameshift += "*";
    } else if (trace == 1) {
        i--;
        j--;
        k = i + 1;
        final_seq1 += translate_codon(str + k - 1);
        final_seq2 += ref_seq[j];
        frameshift += " ";
    }
}

void traceV2_1d(char* input_seq, char* ref_seq, int* sc_mat, int* t_sc_mat, int N, int M, int index,
                int* indexes) {
    int i_max = 0, j_max = 0, i = 0, j = 0, max_score = 0, curr_score, print = 0;
    string f1, f2, f3, seq_dna, seq_prot, frameshift;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            curr_score = sc_mat[i * M + j];
            if (curr_score > max_score) {
                max_score = curr_score;
                i_max = i;
                j_max = j;
            }
        }
    }

    i = i_max;
    j = j_max;

    indexes[0] = max_score;
    indexes[1] = i_max;
    indexes[2] = j_max;

    while (sc_mat[i * M + j] != 0) {
        indexes[3] = j;
        routine(t_sc_mat[i * M + j], i, j, input_seq, ref_seq, seq_dna, seq_prot, frameshift);
        print = i;
    }
    cout << endl;
    reverse(seq_dna.begin(), seq_dna.end());
    reverse(seq_prot.begin(), seq_prot.end());
    reverse(frameshift.begin(), frameshift.end());

    three_frame(input_seq, &f1, &f2, &f3);

    cout << "DNA Input:\t";
    for (i = 0; i < strlen(input_seq); i++) {
        cout << input_seq[i];
    }
    cout << endl;
    cout << "\t\t ";
    for (i = 0; i < f1.length(); i++) {
        cout << f1[i] << "  ";
    }
    cout << endl;
    cout << "\t\t  ";
    for (i = 0; i < f2.length(); i++) {
        cout << f2[i] << "  ";
    }
    cout << endl;
    cout << "\t\t   ";
    for (i = 0; i < f3.length(); i++) {
        cout << f3[i] << "  ";
    }
    cout << endl << endl;

    cout << "Frame Match: \t ";
    for (i = 0; i < seq_dna.length(); i++) {
        cout << seq_dna[i] << "  ";
    }
    cout << endl;
    cout << "Frameshift: \t ";
    for (i = 0; i < frameshift.length(); i++) {
        cout << frameshift[i] << "  ";
    }
    cout << endl;
    cout << "Output Prot: \t ";
    for (i = 0; i < seq_prot.length(); i++) {
        cout << seq_prot[i] << "  ";
    }
    cout << endl;
}

void traceV2_1d_check(char* input_seq, char* ref_seq, int* sc_mat, int* t_sc_mat, size_t N,
                      size_t M, int index, int* indexes) {
    N = (int)N;
    M = (int)M;
    int i_max = 0, j_max = 0, max_score = 0, curr_score, i = 0, j = 0;
    string seq_dna, seq_prot, frameshift;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            curr_score = sc_mat[i * M + j];
            if (curr_score > max_score) {
                max_score = curr_score;
                i_max = i;
                j_max = j;
            }
        }
    }

    indexes[0] = max_score;
    indexes[1] = i_max;
    indexes[2] = j_max;

    i = i_max;
    j = j_max;

    myArray[index][0] = to_string(max_score);
    myArray[index][2] = to_string(j_max);

    while (sc_mat[i * M + j] != 0) {
        myArray[index][1] = to_string(j);
        indexes[3] = j;
        routine(t_sc_mat[i * M + j], i, j, input_seq, ref_seq, seq_dna, seq_prot, frameshift);
    }
}

void write_to_excel(int n, int i) {
    string filename = "outputRun6F" + std::to_string(n) + ".csv";
    std::ofstream file(filename, std::ios::app);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing." << std::endl;
        return;
    }
    file << i << "," << myArray[i][4] << "," << myArray[i][0] << "," << myArray[i][1] << ","
         << myArray[i][2] << "," << myArray[i][3] << "\n";

    file.close();
    cout << "Data " << i << " successfully written to output.csv" << endl;
}

int main() {
    int mode, frame, top, device;
    vector<string> dnaInputs, proteinInputs, proteinIdInputs;
    string protein_sequence, DNA_sequence, DNA_sequence_r, file_name;
    GpuTimer timer;

    checkCudaErrors(cudaGetDevice(&device));
    checkCudaErrors(cudaSetDevice(device));

    readBlosum62();
    cudaMemcpyToSymbol(d_blosum62mat, blosum62mat, sizeof(blosum62mat), 0, cudaMemcpyHostToDevice);
    initialize_protein_to_int_map();

    do {
        file_name.clear();
        cout << "Input file name for protein(exclude .fasta/.fastq):" << endl << "--> ";
        cin >> file_name;

        proteinInputs = readFastaSequences(file_name);
        proteinIdInputs = readFastaIDs(file_name);
        if (!proteinInputs.empty() && !proteinIdInputs.empty()) {
            break;
        }
    } while (true);

    do {
        file_name.clear();
        cout << "Input file name for DNA (exclude .fasta/.fastq):" << endl << "--> ";
        cin >> file_name;

        dnaInputs = readFastaSequences(file_name);
        if (!dnaInputs.empty()) {
            break;
        }
    } while (true);

    do {
        cout << "Choose mode [0 - Sequential] [1 - CUDA]:" << endl << "--> ";
        cin >> mode;

        if (mode == 0 || mode == 1)
            break;

        cout << "Invalid input please try again." << endl;
    } while (true);

    do {
        cout << "Choose frame count [3 - 3 Frame] [6 - 6 Frame]:" << endl << "--> ";
        cin >> frame;

        if (frame == 3 || frame == 6)
            break;

        cout << "Invalid input please try again." << endl;
    } while (true);

    do {
        cout << "Show Score [0 - Top 1] [1 - Top 5]:" << endl << "--> ";
        cin >> top;

        if (top == 0 || top == 1)
            break;

        cout << "Invalid input please try again." << endl;

    } while (true);

    double total_r = 0.0;

    for (int index_dna = 0; index_dna < dnaInputs.size(); index_dna++) {

        int* top_scores = new int[5]{};
        int* top_i = new int[5]{};
        int* top_j = new int[5]{};
        int* top_indexes = new int[5]{};
        int* index = new int[4]{};
        int* index_r = new int[4]{};
        int top_hold = 0;

        steady_clock::time_point start, end;

        const int num_prots = static_cast<int>(min<size_t>(proteinInputs.size(), 4096));
        for (int index_prot = 0; index_prot < num_prots; index_prot++) {

            DNA_sequence = dnaInputs[index_dna];
            protein_sequence = proteinInputs[index_prot];
            DNA_sequence_r = reverse_complement(DNA_sequence);

            size_t N = DNA_sequence.length();
            size_t M = protein_sequence.length() + 1;

            size_t N_size = (N) * sizeof(char);
            size_t M_size = (M) * sizeof(char);
            size_t size = (N) * (M) * sizeof(int);

            myArray[index_prot][4] = to_string(protein_sequence.length());

            int* sc_mat = (int*)malloc(size);
            int* ins_mat = (int*)malloc(size);
            int* del_mat = (int*)malloc(size);
            int* sc_mat_hold = (int*)malloc(size);

            int* t_sc_mat = (int*)malloc(size);
            int* t_ins_mat = (int*)malloc(size);
            int* t_del_mat = (int*)malloc(size);
            int* t_sc_mat_hold = (int*)malloc(size);

            int* sc_mat_r = (int*)malloc(size);
            int* ins_mat_r = (int*)malloc(size);
            int* del_mat_r = (int*)malloc(size);
            int* sc_mat_hold_r = (int*)malloc(size);

            int* t_sc_mat_r = (int*)malloc(size);
            int* t_ins_mat_r = (int*)malloc(size);
            int* t_del_mat_r = (int*)malloc(size);
            int* t_sc_mat_hold_r = (int*)malloc(size);

            char* c_DNA_sequence = new char[N_size];
            char* c_protein_sequence = new char[M_size];
            char* c_DNA_sequence_r = new char[N_size];

            memcpy(c_DNA_sequence, DNA_sequence.c_str(), N_size);
            memcpy(c_protein_sequence, protein_sequence.c_str(), M_size);
            memcpy(c_DNA_sequence_r, DNA_sequence_r.c_str(), N_size);

            init_local_v2(c_DNA_sequence, c_protein_sequence, sc_mat, ins_mat, del_mat, t_sc_mat,
                          t_ins_mat, t_del_mat, N, M);
            if (frame == 6)
                init_local_v2(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, ins_mat_r, del_mat_r,
                              t_sc_mat_r, t_ins_mat_r, t_del_mat_r, N, M);

            if (mode == 0) {
                start = steady_clock::now();
                scoring_local_v2(c_DNA_sequence, c_protein_sequence, sc_mat, ins_mat, del_mat,
                                 t_sc_mat, t_ins_mat, t_del_mat, N, M);
                end = steady_clock::now();
                double elapsed1 = duration<double, milli>(end - start).count();

                traceV2_1d_check(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                 index_prot, index);
                myArray[index_prot][3] = to_string(elapsed1);
                top5(index[0], index_prot, index[1], index[2], top_scores, top_i, top_j,
                     top_indexes);
                // write_to_excel(index_dna, index_prot);
                total_r += elapsed1;
                cout << "Run DNA: " << index_dna << " Prot: " << index_prot << endl
                     << "Time in ms: " << elapsed1 << endl;
                cout << "Total runtime: " << total_r << endl;

                if (frame == 6) {
                    start = steady_clock::now();
                    scoring_local_v2(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, ins_mat_r,
                                     del_mat_r, t_sc_mat_r, t_ins_mat_r, t_del_mat_r, N, M);
                    end = steady_clock::now();
                    double elapsed = duration<double, milli>(end - start).count();

                    traceV2_1d_check(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, t_sc_mat_r, N,
                                     M, index_prot, index_r);
                    myArray[index_prot][3] = to_string(elapsed + elapsed1);
                    cout << "Run DNA: " << index_dna << " Prot: " << index_prot << endl
                         << "Time in ms: " << elapsed + elapsed1 << endl;
                    total_r += elapsed;
                    cout << "Total runtime: " << total_r << endl;

                    if (index[0] >= index_r[0]) {
                        myArray[index_prot][0] = to_string(index[0]);
                        myArray[index_prot][1] = to_string(index[3]);
                        myArray[index_prot][2] = to_string(index[2]);
                    } else {
                        myArray[index_prot][0] = to_string(index_r[0]);
                        myArray[index_prot][1] = to_string(index_r[3]);
                        myArray[index_prot][2] = to_string(index_r[2]);
                        top5(index_r[0], index_prot, index_r[1], index_r[2], top_scores, top_i,
                             top_j, top_indexes);
                    }
                    // write_to_excel(index_dna, index_prot);
                }
            } else if (mode == 1) {

                char* d_DNA_sequence;
                char* d_protein_sequence;
                char* d_DNA_sequence_r;

                int* d_sc_mat;
                int* d_ins_mat;
                int* d_del_mat;

                int* d_t_sc_mat;
                int* d_t_ins_mat;
                int* d_t_del_mat;

                int* d_sc_mat_r;
                int* d_ins_mat_r;
                int* d_del_mat_r;

                int* d_t_sc_mat_r;
                int* d_t_ins_mat_r;
                int* d_t_del_mat_r;

                cudaStream_t stream1, stream2;
                checkCudaErrors(cudaStreamCreate(&stream1));
                if (frame == 6)
                    checkCudaErrors(cudaStreamCreate(&stream2));

                checkCudaErrors(cudaMalloc(&d_DNA_sequence, N_size));
                checkCudaErrors(cudaMalloc(&d_protein_sequence, M_size));
                checkCudaErrors(cudaMalloc(&d_DNA_sequence_r, N_size));

                checkCudaErrors(
                    cudaMemcpy(d_DNA_sequence, c_DNA_sequence, N_size, cudaMemcpyHostToDevice));
                checkCudaErrors(cudaMemcpy(d_protein_sequence, c_protein_sequence, M_size,
                                           cudaMemcpyHostToDevice));
                checkCudaErrors(
                    cudaMemcpy(d_DNA_sequence_r, c_DNA_sequence_r, N_size, cudaMemcpyHostToDevice));

                checkCudaErrors(cudaMalloc(&d_sc_mat, size));
                checkCudaErrors(cudaMalloc(&d_ins_mat, size));
                checkCudaErrors(cudaMalloc(&d_del_mat, size));

                checkCudaErrors(cudaMalloc(&d_t_sc_mat, size));
                checkCudaErrors(cudaMalloc(&d_t_ins_mat, size));
                checkCudaErrors(cudaMalloc(&d_t_del_mat, size));

                checkCudaErrors(cudaMalloc(&d_sc_mat_r, size));
                checkCudaErrors(cudaMalloc(&d_ins_mat_r, size));
                checkCudaErrors(cudaMalloc(&d_del_mat_r, size));

                checkCudaErrors(cudaMalloc(&d_t_sc_mat_r, size));
                checkCudaErrors(cudaMalloc(&d_t_ins_mat_r, size));
                checkCudaErrors(cudaMalloc(&d_t_del_mat_r, size));

                checkCudaErrors(cudaMemset(d_sc_mat, 0, size));
                checkCudaErrors(cudaMemset(d_ins_mat, 0, size));
                checkCudaErrors(cudaMemset(d_del_mat, 0, size));
                checkCudaErrors(cudaMemset(d_t_sc_mat, 0, size));
                checkCudaErrors(cudaMemset(d_t_ins_mat, 0, size));
                checkCudaErrors(cudaMemset(d_t_del_mat, 0, size));

                checkCudaErrors(cudaMemset(d_sc_mat_r, 0, size));
                checkCudaErrors(cudaMemset(d_ins_mat_r, 0, size));
                checkCudaErrors(cudaMemset(d_del_mat_r, 0, size));
                checkCudaErrors(cudaMemset(d_t_sc_mat_r, 0, size));
                checkCudaErrors(cudaMemset(d_t_ins_mat_r, 0, size));
                checkCudaErrors(cudaMemset(d_t_del_mat_r, 0, size));

                checkCudaErrors(cudaMemcpy(d_sc_mat, sc_mat, size, cudaMemcpyHostToDevice));
                checkCudaErrors(cudaMemcpy(d_ins_mat, ins_mat, size, cudaMemcpyHostToDevice));
                checkCudaErrors(cudaMemcpy(d_del_mat, del_mat, size, cudaMemcpyHostToDevice));
                checkCudaErrors(cudaMemcpy(d_t_sc_mat, t_sc_mat, size, cudaMemcpyHostToDevice));
                checkCudaErrors(cudaMemcpy(d_t_ins_mat, t_ins_mat, size, cudaMemcpyHostToDevice));
                checkCudaErrors(cudaMemcpy(d_t_del_mat, t_del_mat, size, cudaMemcpyHostToDevice));

                checkCudaErrors(cudaMemcpy(d_sc_mat_r, sc_mat_r, size, cudaMemcpyHostToDevice));
                checkCudaErrors(cudaMemcpy(d_ins_mat_r, ins_mat_r, size, cudaMemcpyHostToDevice));
                checkCudaErrors(cudaMemcpy(d_del_mat_r, del_mat_r, size, cudaMemcpyHostToDevice));
                checkCudaErrors(cudaMemcpy(d_t_sc_mat_r, t_sc_mat_r, size, cudaMemcpyHostToDevice));
                checkCudaErrors(
                    cudaMemcpy(d_t_ins_mat_r, t_ins_mat_r, size, cudaMemcpyHostToDevice));
                checkCudaErrors(
                    cudaMemcpy(d_t_del_mat_r, t_del_mat_r, size, cudaMemcpyHostToDevice));

                dim3 blockDimMain(32, 32);

                unsigned int submatrixSide = blockDimMain.y;
                unsigned int numSubmatrixRows =
                    ((unsigned int)N + submatrixSide - 1) / submatrixSide;
                unsigned int numSubmatrixCols =
                    ((unsigned int)M + submatrixSide - 1) / submatrixSide;

                if (frame == 3) {

                    timer.Start();

                    for (unsigned int diag = 0; diag < numSubmatrixRows + numSubmatrixCols - 1;
                         ++diag) {
                        unsigned int first_sub_y_curr_diag =
                            std::max(0, (int)diag - (int)(numSubmatrixCols - 1));
                        unsigned int last_sub_y_curr_diag =
                            std::min((int)diag, (int)numSubmatrixRows - 1);

                        unsigned int curr_diag_blocks = 0;
                        if (last_sub_y_curr_diag >= first_sub_y_curr_diag)
                            curr_diag_blocks = last_sub_y_curr_diag - first_sub_y_curr_diag + 1;
                        dim3 currGridDim(1, 1);

                        if (curr_diag_blocks > 0) {
                            currGridDim = dim3(curr_diag_blocks, 1);
                        }
                        scoring_local_v2_cuda<<<currGridDim, blockDimMain, 0, stream1>>>(
                            d_DNA_sequence, d_protein_sequence, d_sc_mat, d_ins_mat, d_del_mat,
                            d_t_sc_mat, d_t_ins_mat, d_t_del_mat, (unsigned int)N, (unsigned int)M,
                            diag, submatrixSide, first_sub_y_curr_diag);
                    }
                    checkCudaErrors(cudaStreamSynchronize(stream1));

                    timer.Stop();

                    checkCudaErrors(cudaStreamDestroy(stream1));
                    checkCudaErrors(cudaMemcpy(sc_mat, d_sc_mat, size, cudaMemcpyDeviceToHost));
                    checkCudaErrors(cudaMemcpy(t_sc_mat, d_t_sc_mat, size, cudaMemcpyDeviceToHost));

                    traceV2_1d_check(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                     index_prot, index);
                    cout << "Run DNA: " << index_dna << " Prot: " << index_prot << endl
                         << "Time in ms: " << timer.Elapsed() << endl;
                    myArray[index_prot][3] = to_string(timer.Elapsed());
                    // write_to_excel(index_dna, index_prot);
                    top5(index[0], index_prot, index[1], index[2], top_scores, top_i, top_j,
                         top_indexes);
                    total_r += timer.Elapsed();
                    cout << "Total runtime: " << total_r << endl;
                }

                if (frame == 6) {

                    timer.Start();

                    for (unsigned int diag = 0; diag < numSubmatrixRows + numSubmatrixCols - 1;
                         ++diag) {
                        unsigned int first_sub_y_curr_diag =
                            std::max(0, (int)diag - (int)(numSubmatrixCols - 1));
                        unsigned int last_sub_y_curr_diag =
                            std::min((int)diag, (int)numSubmatrixRows - 1);

                        unsigned int curr_diag_blocks = 0;
                        if (last_sub_y_curr_diag >= first_sub_y_curr_diag)
                            curr_diag_blocks = last_sub_y_curr_diag - first_sub_y_curr_diag + 1;
                        dim3 currGridDim(1, 1);

                        if (curr_diag_blocks > 0) {
                            currGridDim = dim3(curr_diag_blocks, 1);
                        }
                        scoring_local_v2_cuda<<<currGridDim, blockDimMain, 0, stream1>>>(
                            d_DNA_sequence, d_protein_sequence, d_sc_mat, d_ins_mat, d_del_mat,
                            d_t_sc_mat, d_t_ins_mat, d_t_del_mat, (unsigned int)N, (unsigned int)M,
                            diag, submatrixSide, first_sub_y_curr_diag);
                        scoring_local_v2_cuda<<<currGridDim, blockDimMain, 0, stream2>>>(
                            d_DNA_sequence_r, d_protein_sequence, d_sc_mat_r, d_ins_mat_r,
                            d_del_mat_r, d_t_sc_mat_r, d_t_ins_mat_r, d_t_del_mat_r,
                            (unsigned int)N, (unsigned int)M, diag, submatrixSide,
                            first_sub_y_curr_diag);
                    }
                    checkCudaErrors(cudaStreamSynchronize(stream1));
                    checkCudaErrors(cudaStreamSynchronize(stream2));

                    timer.Stop();

                    checkCudaErrors(cudaStreamDestroy(stream1));
                    checkCudaErrors(cudaStreamDestroy(stream2));

                    checkCudaErrors(cudaMemcpy(sc_mat, d_sc_mat, size, cudaMemcpyDeviceToHost));
                    checkCudaErrors(cudaMemcpy(t_sc_mat, d_t_sc_mat, size, cudaMemcpyDeviceToHost));

                    checkCudaErrors(cudaMemcpy(sc_mat_r, d_sc_mat_r, size, cudaMemcpyDeviceToHost));
                    checkCudaErrors(
                        cudaMemcpy(t_sc_mat_r, d_t_sc_mat_r, size, cudaMemcpyDeviceToHost));

                    traceV2_1d_check(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                     index_prot, index);
                    cout << "Run DNA: " << index_dna << " Prot: " << index_prot << endl
                         << "Time in ms: " << timer.Elapsed() << endl;
                    myArray[index_prot][3] = to_string(timer.Elapsed());

                    traceV2_1d_check(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, t_sc_mat_r, N,
                                     M, index_prot, index_r);
                    cout << "Run DNA: " << index_dna << " Prot: " << index_prot << endl
                         << "Time in ms: " << timer.Elapsed() << endl;

                    total_r += timer.Elapsed();
                    cout << "Total runtime: " << total_r << endl;

                    if (index[0] >= index_r[0]) {
                        myArray[index_prot][0] = to_string(index[0]);
                        myArray[index_prot][1] = to_string(index[3]);
                        myArray[index_prot][2] = to_string(index[2]);
                        top5(index[0], index_prot, index[1], index[2], top_scores, top_i, top_j,
                             top_indexes);
                    } else {
                        myArray[index_prot][0] = to_string(index_r[0]);
                        myArray[index_prot][1] = to_string(index_r[3]);
                        myArray[index_prot][2] = to_string(index_r[2]);
                        top5(index_r[0], index_prot, index_r[1], index_r[2], top_scores, top_i,
                             top_j, top_indexes);
                    }

                    // write_to_excel(index_dna, index_prot);
                }

                checkCudaErrors(cudaFree(d_DNA_sequence));
                checkCudaErrors(cudaFree(d_protein_sequence));
                checkCudaErrors(cudaFree(d_DNA_sequence_r));

                checkCudaErrors(cudaFree(d_sc_mat));
                checkCudaErrors(cudaFree(d_ins_mat));
                checkCudaErrors(cudaFree(d_del_mat));

                checkCudaErrors(cudaFree(d_t_sc_mat));
                checkCudaErrors(cudaFree(d_t_ins_mat));
                checkCudaErrors(cudaFree(d_t_del_mat));

                checkCudaErrors(cudaFree(d_sc_mat_r));
                checkCudaErrors(cudaFree(d_ins_mat_r));
                checkCudaErrors(cudaFree(d_del_mat_r));

                checkCudaErrors(cudaFree(d_t_sc_mat_r));
                checkCudaErrors(cudaFree(d_t_ins_mat_r));
                checkCudaErrors(cudaFree(d_t_del_mat_r));
                cudaDeviceSynchronize();
            }

            free(sc_mat);
            free(ins_mat);
            free(del_mat);
            free(sc_mat_hold);

            free(t_sc_mat);
            free(t_ins_mat);
            free(t_del_mat);
            free(t_sc_mat_hold);

            free(sc_mat_r);
            free(ins_mat_r);
            free(del_mat_r);
            free(sc_mat_hold_r);

            free(t_sc_mat_r);
            free(t_ins_mat_r);
            free(t_del_mat_r);
            free(t_sc_mat_hold_r);

            delete[] c_protein_sequence;
            delete[] c_DNA_sequence;
            delete[] c_DNA_sequence_r;
        }

        if (top == 0) {
            for (int i = 0; i < 5; i++) {
                DNA_sequence = dnaInputs[0];
                protein_sequence = proteinInputs[top_indexes[i]];
                DNA_sequence_r = reverse_complement(DNA_sequence);

                size_t N = DNA_sequence.length();
                size_t M = protein_sequence.length() + 1;

                size_t N_size = (N) * sizeof(char);
                size_t M_size = (M) * sizeof(char);
                size_t size = (N) * (M) * sizeof(int);
                top_hold = top_scores[0];

                int* sc_mat = (int*)malloc(size);
                int* ins_mat = (int*)malloc(size);
                int* del_mat = (int*)malloc(size);
                int* sc_mat_hold = (int*)malloc(size);

                int* t_sc_mat = (int*)malloc(size);
                int* t_ins_mat = (int*)malloc(size);
                int* t_del_mat = (int*)malloc(size);
                int* t_sc_mat_hold = (int*)malloc(size);

                int* sc_mat_r = (int*)malloc(size);
                int* ins_mat_r = (int*)malloc(size);
                int* del_mat_r = (int*)malloc(size);
                int* sc_mat_hold_r = (int*)malloc(size);

                int* t_sc_mat_r = (int*)malloc(size);
                int* t_ins_mat_r = (int*)malloc(size);
                int* t_del_mat_r = (int*)malloc(size);
                int* t_sc_mat_hold_r = (int*)malloc(size);

                char* c_DNA_sequence = new char[N_size];
                char* c_protein_sequence = new char[M_size];
                char* c_DNA_sequence_r = new char[N_size];

                memcpy(c_DNA_sequence, DNA_sequence.c_str(), N_size);
                memcpy(c_protein_sequence, protein_sequence.c_str(), M_size);
                memcpy(c_DNA_sequence_r, DNA_sequence_r.c_str(), N_size);

                init_local_v2(c_DNA_sequence, c_protein_sequence, sc_mat, ins_mat, del_mat,
                              t_sc_mat, t_ins_mat, t_del_mat, N, M);
                if (frame == 6)
                    init_local_v2(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, ins_mat_r,
                                  del_mat_r, t_sc_mat_r, t_ins_mat_r, t_del_mat_r, N, M);

                if (mode == 0 && top_hold == top_scores[i]) {
                    cout << proteinIdInputs[top_indexes[i]] << endl
                         << "Sequence: " << top_indexes[i] << endl;
                    start = steady_clock::now();
                    scoring_local_v2(c_DNA_sequence, c_protein_sequence, sc_mat, ins_mat, del_mat,
                                     t_sc_mat, t_ins_mat, t_del_mat, N, M);
                    end = steady_clock::now();
                    double elapsed1 = duration<double, milli>(end - start).count();

                    traceV2_1d_check(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                     top_indexes[i], index);

                    if (frame == 3) {
                        traceV2_1d(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                   top_indexes[i], index);
                        cout << endl << "Score: " << top_scores[i] << endl;
                        cout << "Start to End match in Protein: " << myArray[top_indexes[i]][1]
                             << "-" << myArray[top_indexes[i]][2] << endl
                             << endl;
                        cout << "Time in ms: " << elapsed1 << endl << endl;
                    }

                    if (frame == 6) {
                        start = steady_clock::now();
                        scoring_local_v2(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, ins_mat_r,
                                         del_mat_r, t_sc_mat_r, t_ins_mat_r, t_del_mat_r, N, M);
                        end = steady_clock::now();
                        double elapsed = duration<double, milli>(end - start).count();

                        traceV2_1d_check(c_DNA_sequence, c_protein_sequence, sc_mat_r, t_sc_mat_r,
                                         N, M, top_indexes[i], index_r);

                        if (index[0] >= index_r[0]) {
                            traceV2_1d(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                       top_indexes[i], index);
                            cout << endl << "Score: " << top_scores[i] << endl;
                            cout << "Start to End match in Protein: " << index[3] << "-" << index[2]
                                 << endl
                                 << endl;
                            cout << "Time in ms: " << elapsed1 << endl << endl;
                        } else {
                            cout << "Reverse: " << endl;
                            traceV2_1d(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, t_sc_mat_r,
                                       N, M, top_indexes[i], index_r);
                            cout << endl << "Score: " << top_scores[i] << endl;
                            cout << "Start to End match in Protein: " << index_r[3] << "-"
                                 << index_r[2] << endl
                                 << endl;
                            cout << "Time in ms: " << elapsed << endl << endl;
                        }
                    }
                } else if (mode == 1 && top_hold == top_scores[i]) {
                    cout << proteinIdInputs[top_indexes[i]] << endl
                         << "Sequence: " << top_indexes[i] << endl;
                    char* d_DNA_sequence;
                    char* d_protein_sequence;
                    char* d_DNA_sequence_r;

                    int* d_sc_mat;
                    int* d_ins_mat;
                    int* d_del_mat;

                    int* d_t_sc_mat;
                    int* d_t_ins_mat;
                    int* d_t_del_mat;

                    int* d_sc_mat_r;
                    int* d_ins_mat_r;
                    int* d_del_mat_r;

                    int* d_t_sc_mat_r;
                    int* d_t_ins_mat_r;
                    int* d_t_del_mat_r;

                    cudaStream_t stream1, stream2;
                    checkCudaErrors(cudaStreamCreate(&stream1));
                    if (frame == 6)
                        checkCudaErrors(cudaStreamCreate(&stream2));

                    checkCudaErrors(cudaMalloc(&d_DNA_sequence, N_size));
                    checkCudaErrors(cudaMalloc(&d_protein_sequence, M_size));
                    checkCudaErrors(cudaMalloc(&d_DNA_sequence_r, N_size));

                    checkCudaErrors(
                        cudaMemcpy(d_DNA_sequence, c_DNA_sequence, N_size, cudaMemcpyHostToDevice));
                    checkCudaErrors(cudaMemcpy(d_protein_sequence, c_protein_sequence, M_size,
                                               cudaMemcpyHostToDevice));
                    checkCudaErrors(cudaMemcpy(d_DNA_sequence_r, c_DNA_sequence_r, N_size,
                                               cudaMemcpyHostToDevice));

                    checkCudaErrors(cudaMalloc(&d_sc_mat, size));
                    checkCudaErrors(cudaMalloc(&d_ins_mat, size));
                    checkCudaErrors(cudaMalloc(&d_del_mat, size));

                    checkCudaErrors(cudaMalloc(&d_t_sc_mat, size));
                    checkCudaErrors(cudaMalloc(&d_t_ins_mat, size));
                    checkCudaErrors(cudaMalloc(&d_t_del_mat, size));

                    checkCudaErrors(cudaMalloc(&d_sc_mat_r, size));
                    checkCudaErrors(cudaMalloc(&d_ins_mat_r, size));
                    checkCudaErrors(cudaMalloc(&d_del_mat_r, size));

                    checkCudaErrors(cudaMalloc(&d_t_sc_mat_r, size));
                    checkCudaErrors(cudaMalloc(&d_t_ins_mat_r, size));
                    checkCudaErrors(cudaMalloc(&d_t_del_mat_r, size));

                    checkCudaErrors(cudaMemset(d_sc_mat, 0, size));
                    checkCudaErrors(cudaMemset(d_ins_mat, 0, size));
                    checkCudaErrors(cudaMemset(d_del_mat, 0, size));
                    checkCudaErrors(cudaMemset(d_t_sc_mat, 0, size));
                    checkCudaErrors(cudaMemset(d_t_ins_mat, 0, size));
                    checkCudaErrors(cudaMemset(d_t_del_mat, 0, size));

                    checkCudaErrors(cudaMemset(d_sc_mat_r, 0, size));
                    checkCudaErrors(cudaMemset(d_ins_mat_r, 0, size));
                    checkCudaErrors(cudaMemset(d_del_mat_r, 0, size));
                    checkCudaErrors(cudaMemset(d_t_sc_mat_r, 0, size));
                    checkCudaErrors(cudaMemset(d_t_ins_mat_r, 0, size));
                    checkCudaErrors(cudaMemset(d_t_del_mat_r, 0, size));

                    checkCudaErrors(cudaMemcpy(d_sc_mat, sc_mat, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(cudaMemcpy(d_ins_mat, ins_mat, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(cudaMemcpy(d_del_mat, del_mat, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(cudaMemcpy(d_t_sc_mat, t_sc_mat, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(
                        cudaMemcpy(d_t_ins_mat, t_ins_mat, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(
                        cudaMemcpy(d_t_del_mat, t_del_mat, size, cudaMemcpyHostToDevice));

                    checkCudaErrors(cudaMemcpy(d_sc_mat_r, sc_mat_r, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(
                        cudaMemcpy(d_ins_mat_r, ins_mat_r, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(
                        cudaMemcpy(d_del_mat_r, del_mat_r, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(
                        cudaMemcpy(d_t_sc_mat_r, t_sc_mat_r, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(
                        cudaMemcpy(d_t_ins_mat_r, t_ins_mat_r, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(
                        cudaMemcpy(d_t_del_mat_r, t_del_mat_r, size, cudaMemcpyHostToDevice));

                    dim3 blockDimMain(32, 32);

                    unsigned int submatrixSide = blockDimMain.y;
                    unsigned int numSubmatrixRows =
                        ((unsigned int)N + submatrixSide - 1) / submatrixSide;
                    unsigned int numSubmatrixCols =
                        ((unsigned int)M + submatrixSide - 1) / submatrixSide;

                    if (frame == 3) {

                        timer.Start();

                        for (unsigned int diag = 0; diag < numSubmatrixRows + numSubmatrixCols - 1;
                             ++diag) {
                            unsigned int first_sub_y_curr_diag =
                                std::max(0, (int)diag - (int)(numSubmatrixCols - 1));
                            unsigned int last_sub_y_curr_diag =
                                std::min((int)diag, (int)numSubmatrixRows - 1);

                            unsigned int curr_diag_blocks = 0;
                            if (last_sub_y_curr_diag >= first_sub_y_curr_diag)
                                curr_diag_blocks = last_sub_y_curr_diag - first_sub_y_curr_diag + 1;
                            dim3 currGridDim(1, 1);

                            if (curr_diag_blocks > 0) {
                                currGridDim = dim3(curr_diag_blocks, 1);
                            }
                            scoring_local_v2_cuda<<<currGridDim, blockDimMain, 0, stream1>>>(
                                d_DNA_sequence, d_protein_sequence, d_sc_mat, d_ins_mat, d_del_mat,
                                d_t_sc_mat, d_t_ins_mat, d_t_del_mat, (unsigned int)N,
                                (unsigned int)M, diag, submatrixSide, first_sub_y_curr_diag);
                        }
                        checkCudaErrors(cudaStreamSynchronize(stream1));

                        timer.Stop();

                        checkCudaErrors(cudaStreamDestroy(stream1));
                        checkCudaErrors(cudaMemcpy(sc_mat, d_sc_mat, size, cudaMemcpyDeviceToHost));
                        checkCudaErrors(
                            cudaMemcpy(t_sc_mat, d_t_sc_mat, size, cudaMemcpyDeviceToHost));

                        traceV2_1d(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                   top_indexes[i], index);
                        cout << endl << "Score: " << top_scores[i] << endl;
                        cout << "Start to End match in Protein: " << myArray[top_indexes[i]][1]
                             << "-" << myArray[top_indexes[i]][2] << endl
                             << endl;
                        cout << "Time in ms: " << timer.Elapsed() << endl << endl;
                    }

                    if (frame == 6) {

                        timer.Start();

                        for (unsigned int diag = 0; diag < numSubmatrixRows + numSubmatrixCols - 1;
                             ++diag) {
                            unsigned int first_sub_y_curr_diag =
                                std::max(0, (int)diag - (int)(numSubmatrixCols - 1));
                            unsigned int last_sub_y_curr_diag =
                                std::min((int)diag, (int)numSubmatrixRows - 1);

                            unsigned int curr_diag_blocks = 0;
                            if (last_sub_y_curr_diag >= first_sub_y_curr_diag)
                                curr_diag_blocks = last_sub_y_curr_diag - first_sub_y_curr_diag + 1;
                            dim3 currGridDim(1, 1);

                            if (curr_diag_blocks > 0) {
                                currGridDim = dim3(curr_diag_blocks, 1);
                            }
                            scoring_local_v2_cuda<<<currGridDim, blockDimMain, 0, stream1>>>(
                                d_DNA_sequence, d_protein_sequence, d_sc_mat, d_ins_mat, d_del_mat,
                                d_t_sc_mat, d_t_ins_mat, d_t_del_mat, (unsigned int)N,
                                (unsigned int)M, diag, submatrixSide, first_sub_y_curr_diag);
                            scoring_local_v2_cuda<<<currGridDim, blockDimMain, 0, stream2>>>(
                                d_DNA_sequence_r, d_protein_sequence, d_sc_mat_r, d_ins_mat_r,
                                d_del_mat_r, d_t_sc_mat_r, d_t_ins_mat_r, d_t_del_mat_r,
                                (unsigned int)N, (unsigned int)M, diag, submatrixSide,
                                first_sub_y_curr_diag);
                        }
                        checkCudaErrors(cudaStreamSynchronize(stream1));
                        checkCudaErrors(cudaStreamSynchronize(stream2));

                        timer.Stop();

                        checkCudaErrors(cudaStreamDestroy(stream1));
                        checkCudaErrors(cudaStreamDestroy(stream2));

                        checkCudaErrors(cudaMemcpy(sc_mat, d_sc_mat, size, cudaMemcpyDeviceToHost));
                        checkCudaErrors(
                            cudaMemcpy(t_sc_mat, d_t_sc_mat, size, cudaMemcpyDeviceToHost));

                        checkCudaErrors(
                            cudaMemcpy(sc_mat_r, d_sc_mat_r, size, cudaMemcpyDeviceToHost));
                        checkCudaErrors(
                            cudaMemcpy(t_sc_mat_r, d_t_sc_mat_r, size, cudaMemcpyDeviceToHost));

                        traceV2_1d_check(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                         top_indexes[i], index);
                        traceV2_1d_check(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, t_sc_mat_r,
                                         N, M, top_indexes[i], index_r);

                        if (index[0] >= index_r[0]) {
                            traceV2_1d(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                       top_indexes[i], index);
                            cout << endl << "Score: " << top_scores[i] << endl;
                            cout << "Start to End match in Protein: " << index[3] << "-" << index[2]
                                 << endl
                                 << endl;
                            cout << "Time in ms: " << timer.Elapsed() << endl << endl;
                        } else {
                            cout << "Reverse: " << endl;
                            traceV2_1d(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, t_sc_mat_r,
                                       N, M, top_indexes[i], index_r);
                            cout << endl << "Score: " << top_scores[i] << endl;
                            cout << "Start to End match in Protein: " << index_r[3] << "-"
                                 << index_r[2] << endl
                                 << endl;
                            cout << "Time in ms: " << timer.Elapsed() << endl << endl;
                        }
                    }

                    checkCudaErrors(cudaFree(d_DNA_sequence));
                    checkCudaErrors(cudaFree(d_protein_sequence));
                    checkCudaErrors(cudaFree(d_DNA_sequence_r));

                    checkCudaErrors(cudaFree(d_sc_mat));
                    checkCudaErrors(cudaFree(d_ins_mat));
                    checkCudaErrors(cudaFree(d_del_mat));

                    checkCudaErrors(cudaFree(d_t_sc_mat));
                    checkCudaErrors(cudaFree(d_t_ins_mat));
                    checkCudaErrors(cudaFree(d_t_del_mat));

                    checkCudaErrors(cudaFree(d_sc_mat_r));
                    checkCudaErrors(cudaFree(d_ins_mat_r));
                    checkCudaErrors(cudaFree(d_del_mat_r));

                    checkCudaErrors(cudaFree(d_t_sc_mat_r));
                    checkCudaErrors(cudaFree(d_t_ins_mat_r));
                    checkCudaErrors(cudaFree(d_t_del_mat_r));
                    cudaDeviceSynchronize();
                }

                free(sc_mat);
                free(ins_mat);
                free(del_mat);
                free(sc_mat_hold);

                free(t_sc_mat);
                free(t_ins_mat);
                free(t_del_mat);
                free(t_sc_mat_hold);

                free(sc_mat_r);
                free(ins_mat_r);
                free(del_mat_r);
                free(sc_mat_hold_r);

                free(t_sc_mat_r);
                free(t_ins_mat_r);
                free(t_del_mat_r);
                free(t_sc_mat_hold_r);

                delete[] c_protein_sequence;
                delete[] c_DNA_sequence;
                delete[] c_DNA_sequence_r;
            }
        } else if (top == 1) {
            for (int i = 0; i < 5; i++) {
                DNA_sequence = dnaInputs[0];
                protein_sequence = proteinInputs[top_indexes[i]];
                DNA_sequence_r = reverse_complement(DNA_sequence);

                size_t N = DNA_sequence.length();
                size_t M = protein_sequence.length() + 1;

                size_t N_size = (N) * sizeof(char);
                size_t M_size = (M) * sizeof(char);
                size_t size = (N) * (M) * sizeof(int);

                int* sc_mat = (int*)malloc(size);
                int* ins_mat = (int*)malloc(size);
                int* del_mat = (int*)malloc(size);
                int* sc_mat_hold = (int*)malloc(size);

                int* t_sc_mat = (int*)malloc(size);
                int* t_ins_mat = (int*)malloc(size);
                int* t_del_mat = (int*)malloc(size);
                int* t_sc_mat_hold = (int*)malloc(size);

                int* sc_mat_r = (int*)malloc(size);
                int* ins_mat_r = (int*)malloc(size);
                int* del_mat_r = (int*)malloc(size);
                int* sc_mat_hold_r = (int*)malloc(size);

                int* t_sc_mat_r = (int*)malloc(size);
                int* t_ins_mat_r = (int*)malloc(size);
                int* t_del_mat_r = (int*)malloc(size);
                int* t_sc_mat_hold_r = (int*)malloc(size);

                char* c_DNA_sequence = new char[N_size];
                char* c_protein_sequence = new char[M_size];
                char* c_DNA_sequence_r = new char[N_size];

                memcpy(c_DNA_sequence, DNA_sequence.c_str(), N_size);
                memcpy(c_protein_sequence, protein_sequence.c_str(), M_size);
                memcpy(c_DNA_sequence_r, DNA_sequence_r.c_str(), N_size);

                init_local_v2(c_DNA_sequence, c_protein_sequence, sc_mat, ins_mat, del_mat,
                              t_sc_mat, t_ins_mat, t_del_mat, N, M);
                if (frame == 6)
                    init_local_v2(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, ins_mat_r,
                                  del_mat_r, t_sc_mat_r, t_ins_mat_r, t_del_mat_r, N, M);

                if (mode == 0) {
                    cout << proteinIdInputs[top_indexes[i]] << endl
                         << "Sequence: " << top_indexes[i] << endl;
                    start = steady_clock::now();
                    scoring_local_v2(c_DNA_sequence, c_protein_sequence, sc_mat, ins_mat, del_mat,
                                     t_sc_mat, t_ins_mat, t_del_mat, N, M);
                    end = steady_clock::now();
                    double elapsed1 = duration<double, milli>(end - start).count();

                    traceV2_1d_check(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                     top_indexes[i], index);

                    if (frame == 3) {
                        traceV2_1d(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                   top_indexes[i], index);
                        cout << endl << "Score: " << top_scores[i] << endl;
                        cout << "Start to End match in Protein: " << index[3] << "-" << index[2]
                             << endl
                             << endl;
                        cout << "Time in ms: " << elapsed1 << endl << endl;
                    }

                    if (frame == 6) {
                        start = steady_clock::now();
                        scoring_local_v2(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, ins_mat_r,
                                         del_mat_r, t_sc_mat_r, t_ins_mat_r, t_del_mat_r, N, M);
                        end = steady_clock::now();
                        double elapsed = duration<double, milli>(end - start).count();

                        traceV2_1d_check(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, t_sc_mat_r,
                                         N, M, top_indexes[i], index_r);

                        if (index[0] >= index_r[0]) {
                            traceV2_1d(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                       top_indexes[i], index);
                            cout << endl << "Score: " << top_scores[i] << endl;
                            cout << "Start to End match in Protein: " << index[3] << "-" << index[2]
                                 << endl
                                 << endl;
                            cout << "Time in ms: " << elapsed1 << endl << endl;
                        } else {
                            cout << "Reverse: " << endl;
                            traceV2_1d(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, t_sc_mat_r,
                                       N, M, top_indexes[i], index_r);
                            cout << endl << "Score: " << top_scores[i] << endl;
                            cout << "Start to End match in Protein: " << index_r[3] << "-"
                                 << index_r[2] << endl
                                 << endl;
                            cout << "Time in ms: " << elapsed << endl << endl;
                        }
                    }
                } else if (mode == 1) {
                    cout << proteinIdInputs[top_indexes[i]] << endl
                         << "Sequence: " << top_indexes[i] << endl;
                    char* d_DNA_sequence;
                    char* d_protein_sequence;
                    char* d_DNA_sequence_r;

                    int* d_sc_mat;
                    int* d_ins_mat;
                    int* d_del_mat;

                    int* d_t_sc_mat;
                    int* d_t_ins_mat;
                    int* d_t_del_mat;

                    int* d_sc_mat_r;
                    int* d_ins_mat_r;
                    int* d_del_mat_r;

                    int* d_t_sc_mat_r;
                    int* d_t_ins_mat_r;
                    int* d_t_del_mat_r;

                    cudaStream_t stream1, stream2;
                    checkCudaErrors(cudaStreamCreate(&stream1));
                    if (frame == 6)
                        checkCudaErrors(cudaStreamCreate(&stream2));

                    checkCudaErrors(cudaMalloc(&d_DNA_sequence, N_size));
                    checkCudaErrors(cudaMalloc(&d_protein_sequence, M_size));
                    checkCudaErrors(cudaMalloc(&d_DNA_sequence_r, N_size));

                    checkCudaErrors(
                        cudaMemcpy(d_DNA_sequence, c_DNA_sequence, N_size, cudaMemcpyHostToDevice));
                    checkCudaErrors(cudaMemcpy(d_protein_sequence, c_protein_sequence, M_size,
                                               cudaMemcpyHostToDevice));
                    checkCudaErrors(cudaMemcpy(d_DNA_sequence_r, c_DNA_sequence_r, N_size,
                                               cudaMemcpyHostToDevice));

                    checkCudaErrors(cudaMalloc(&d_sc_mat, size));
                    checkCudaErrors(cudaMalloc(&d_ins_mat, size));
                    checkCudaErrors(cudaMalloc(&d_del_mat, size));

                    checkCudaErrors(cudaMalloc(&d_t_sc_mat, size));
                    checkCudaErrors(cudaMalloc(&d_t_ins_mat, size));
                    checkCudaErrors(cudaMalloc(&d_t_del_mat, size));

                    checkCudaErrors(cudaMalloc(&d_sc_mat_r, size));
                    checkCudaErrors(cudaMalloc(&d_ins_mat_r, size));
                    checkCudaErrors(cudaMalloc(&d_del_mat_r, size));

                    checkCudaErrors(cudaMalloc(&d_t_sc_mat_r, size));
                    checkCudaErrors(cudaMalloc(&d_t_ins_mat_r, size));
                    checkCudaErrors(cudaMalloc(&d_t_del_mat_r, size));

                    checkCudaErrors(cudaMemset(d_sc_mat, 0, size));
                    checkCudaErrors(cudaMemset(d_ins_mat, 0, size));
                    checkCudaErrors(cudaMemset(d_del_mat, 0, size));
                    checkCudaErrors(cudaMemset(d_t_sc_mat, 0, size));
                    checkCudaErrors(cudaMemset(d_t_ins_mat, 0, size));
                    checkCudaErrors(cudaMemset(d_t_del_mat, 0, size));

                    checkCudaErrors(cudaMemset(d_sc_mat_r, 0, size));
                    checkCudaErrors(cudaMemset(d_ins_mat_r, 0, size));
                    checkCudaErrors(cudaMemset(d_del_mat_r, 0, size));
                    checkCudaErrors(cudaMemset(d_t_sc_mat_r, 0, size));
                    checkCudaErrors(cudaMemset(d_t_ins_mat_r, 0, size));
                    checkCudaErrors(cudaMemset(d_t_del_mat_r, 0, size));

                    checkCudaErrors(cudaMemcpy(d_sc_mat, sc_mat, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(cudaMemcpy(d_ins_mat, ins_mat, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(cudaMemcpy(d_del_mat, del_mat, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(cudaMemcpy(d_t_sc_mat, t_sc_mat, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(
                        cudaMemcpy(d_t_ins_mat, t_ins_mat, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(
                        cudaMemcpy(d_t_del_mat, t_del_mat, size, cudaMemcpyHostToDevice));

                    checkCudaErrors(cudaMemcpy(d_sc_mat_r, sc_mat_r, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(
                        cudaMemcpy(d_ins_mat_r, ins_mat_r, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(
                        cudaMemcpy(d_del_mat_r, del_mat_r, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(
                        cudaMemcpy(d_t_sc_mat_r, t_sc_mat_r, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(
                        cudaMemcpy(d_t_ins_mat_r, t_ins_mat_r, size, cudaMemcpyHostToDevice));
                    checkCudaErrors(
                        cudaMemcpy(d_t_del_mat_r, t_del_mat_r, size, cudaMemcpyHostToDevice));

                    dim3 blockDimMain(32, 32);

                    unsigned int submatrixSide = blockDimMain.y;
                    unsigned int numSubmatrixRows =
                        ((unsigned int)N + submatrixSide - 1) / submatrixSide;
                    unsigned int numSubmatrixCols =
                        ((unsigned int)M + submatrixSide - 1) / submatrixSide;

                    if (frame == 3) {

                        timer.Start();

                        for (unsigned int diag = 0; diag < numSubmatrixRows + numSubmatrixCols - 1;
                             ++diag) {
                            unsigned int first_sub_y_curr_diag =
                                std::max(0, (int)diag - (int)(numSubmatrixCols - 1));
                            unsigned int last_sub_y_curr_diag =
                                std::min((int)diag, (int)numSubmatrixRows - 1);

                            unsigned int curr_diag_blocks = 0;
                            if (last_sub_y_curr_diag >= first_sub_y_curr_diag)
                                curr_diag_blocks = last_sub_y_curr_diag - first_sub_y_curr_diag + 1;
                            dim3 currGridDim(1, 1);

                            if (curr_diag_blocks > 0) {
                                currGridDim = dim3(curr_diag_blocks, 1);
                            }
                            scoring_local_v2_cuda<<<currGridDim, blockDimMain, 0, stream1>>>(
                                d_DNA_sequence, d_protein_sequence, d_sc_mat, d_ins_mat, d_del_mat,
                                d_t_sc_mat, d_t_ins_mat, d_t_del_mat, (unsigned int)N,
                                (unsigned int)M, diag, submatrixSide, first_sub_y_curr_diag);
                        }
                        checkCudaErrors(cudaStreamSynchronize(stream1));

                        timer.Stop();

                        checkCudaErrors(cudaStreamDestroy(stream1));
                        checkCudaErrors(cudaMemcpy(sc_mat, d_sc_mat, size, cudaMemcpyDeviceToHost));
                        checkCudaErrors(
                            cudaMemcpy(t_sc_mat, d_t_sc_mat, size, cudaMemcpyDeviceToHost));

                        traceV2_1d(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                   top_indexes[i], index);
                        cout << endl << "Score: " << top_scores[i] << endl;
                        cout << "Start to End match in Protein: " << index[3] << "-" << index[2]
                             << endl
                             << endl;
                        cout << "Time in ms: " << timer.Elapsed() << endl << endl;
                    }

                    if (frame == 6) {

                        timer.Start();

                        for (unsigned int diag = 0; diag < numSubmatrixRows + numSubmatrixCols - 1;
                             ++diag) {
                            unsigned int first_sub_y_curr_diag =
                                std::max(0, (int)diag - (int)(numSubmatrixCols - 1));
                            unsigned int last_sub_y_curr_diag =
                                std::min((int)diag, (int)numSubmatrixRows - 1);

                            unsigned int curr_diag_blocks = 0;
                            if (last_sub_y_curr_diag >= first_sub_y_curr_diag)
                                curr_diag_blocks = last_sub_y_curr_diag - first_sub_y_curr_diag + 1;
                            dim3 currGridDim(1, 1);

                            if (curr_diag_blocks > 0) {
                                currGridDim = dim3(curr_diag_blocks, 1);
                            }
                            scoring_local_v2_cuda<<<currGridDim, blockDimMain, 0, stream1>>>(
                                d_DNA_sequence, d_protein_sequence, d_sc_mat, d_ins_mat, d_del_mat,
                                d_t_sc_mat, d_t_ins_mat, d_t_del_mat, (unsigned int)N,
                                (unsigned int)M, diag, submatrixSide, first_sub_y_curr_diag);
                            scoring_local_v2_cuda<<<currGridDim, blockDimMain, 0, stream2>>>(
                                d_DNA_sequence_r, d_protein_sequence, d_sc_mat_r, d_ins_mat_r,
                                d_del_mat_r, d_t_sc_mat_r, d_t_ins_mat_r, d_t_del_mat_r,
                                (unsigned int)N, (unsigned int)M, diag, submatrixSide,
                                first_sub_y_curr_diag);
                        }
                        checkCudaErrors(cudaStreamSynchronize(stream1));
                        checkCudaErrors(cudaStreamSynchronize(stream2));

                        timer.Stop();

                        checkCudaErrors(cudaStreamDestroy(stream1));
                        checkCudaErrors(cudaStreamDestroy(stream2));

                        checkCudaErrors(cudaMemcpy(sc_mat, d_sc_mat, size, cudaMemcpyDeviceToHost));
                        checkCudaErrors(
                            cudaMemcpy(t_sc_mat, d_t_sc_mat, size, cudaMemcpyDeviceToHost));

                        checkCudaErrors(
                            cudaMemcpy(sc_mat_r, d_sc_mat_r, size, cudaMemcpyDeviceToHost));
                        checkCudaErrors(
                            cudaMemcpy(t_sc_mat_r, d_t_sc_mat_r, size, cudaMemcpyDeviceToHost));

                        traceV2_1d_check(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                         top_indexes[i], index);
                        traceV2_1d_check(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, t_sc_mat_r,
                                         N, M, top_indexes[i], index_r);

                        if (index[0] >= index_r[0]) {
                            traceV2_1d(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                       top_indexes[i], index);
                            cout << endl << "Score: " << top_scores[i] << endl;
                            cout << "Start to End match in Protein: " << index[3] << "-" << index[2]
                                 << endl
                                 << endl;
                            cout << "Time in ms: " << timer.Elapsed() << endl << endl;
                        } else {
                            cout << "Reverse: " << endl;
                            traceV2_1d(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, t_sc_mat_r,
                                       N, M, top_indexes[i], index_r);
                            cout << endl << "Score: " << top_scores[i] << endl;
                            cout << "Start to End match in Protein: " << index_r[3] << "-"
                                 << index_r[2] << endl
                                 << endl;
                            cout << "Time in ms: " << timer.Elapsed() << endl << endl;
                        }
                    }

                    checkCudaErrors(cudaFree(d_DNA_sequence));
                    checkCudaErrors(cudaFree(d_protein_sequence));
                    checkCudaErrors(cudaFree(d_DNA_sequence_r));

                    checkCudaErrors(cudaFree(d_sc_mat));
                    checkCudaErrors(cudaFree(d_ins_mat));
                    checkCudaErrors(cudaFree(d_del_mat));

                    checkCudaErrors(cudaFree(d_t_sc_mat));
                    checkCudaErrors(cudaFree(d_t_ins_mat));
                    checkCudaErrors(cudaFree(d_t_del_mat));

                    checkCudaErrors(cudaFree(d_sc_mat_r));
                    checkCudaErrors(cudaFree(d_ins_mat_r));
                    checkCudaErrors(cudaFree(d_del_mat_r));

                    checkCudaErrors(cudaFree(d_t_sc_mat_r));
                    checkCudaErrors(cudaFree(d_t_ins_mat_r));
                    checkCudaErrors(cudaFree(d_t_del_mat_r));
                    cudaDeviceSynchronize();
                }

                free(sc_mat);
                free(ins_mat);
                free(del_mat);
                free(sc_mat_hold);

                free(t_sc_mat);
                free(t_ins_mat);
                free(t_del_mat);
                free(t_sc_mat_hold);

                free(sc_mat_r);
                free(ins_mat_r);
                free(del_mat_r);
                free(sc_mat_hold_r);

                free(t_sc_mat_r);
                free(t_ins_mat_r);
                free(t_del_mat_r);
                free(t_sc_mat_hold_r);

                delete[] c_protein_sequence;
                delete[] c_DNA_sequence;
                delete[] c_DNA_sequence_r;
            }
        }

        delete[] top_scores;
        delete[] top_i;
        delete[] top_j;
        delete[] top_indexes;
        delete[] index;
        delete[] index_r;
    }

    cudaFree(d_prot_to_idx);
    cudaFree(d_blosum62mat);

    return 0;
}