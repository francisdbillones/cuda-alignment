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

#include "blosum62.h"
#include "codon.cuh"
#include "cuda_check.cuh"
#include "fasta.h"
#include "gputimer.h"
#include "scoring.cuh"

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

string myArray[40000][5];

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

// An alignment as traceback builds it, from its end back to its start.
struct Alignment {
    string dna, prot, frameshift;
};

void traceV2_1d_check(char* input_seq, char* ref_seq, int* sc_mat, int* t_sc_mat, size_t N,
                      size_t M, int index, int* indexes, Alignment* aln = nullptr) {
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

    if (aln)
        *aln = {seq_dna, seq_prot, frameshift};
}

void traceV2_1d(char* input_seq, char* ref_seq, int* sc_mat, int* t_sc_mat, int N, int M, int index,
                int* indexes) {
    Alignment aln;
    traceV2_1d_check(input_seq, ref_seq, sc_mat, t_sc_mat, N, M, index, indexes, &aln);
    string f1, f2, f3;
    size_t i;

    cout << endl;
    reverse(aln.dna.begin(), aln.dna.end());
    reverse(aln.prot.begin(), aln.prot.end());
    reverse(aln.frameshift.begin(), aln.frameshift.end());

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
    for (i = 0; i < aln.dna.length(); i++) {
        cout << aln.dna[i] << "  ";
    }
    cout << endl;
    cout << "Frameshift: \t ";
    for (i = 0; i < aln.frameshift.length(); i++) {
        cout << aln.frameshift[i] << "  ";
    }
    cout << endl;
    cout << "Output Prot: \t ";
    for (i = 0; i < aln.prot.length(); i++) {
        cout << aln.prot[i] << "  ";
    }
    cout << endl;
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

    read_blosum62();
    initialize_protein_to_int_map();
    upload_blosum62_tables();

    do {
        file_name.clear();
        cout << "Input file name for protein(exclude .fasta/.fastq):" << endl << "--> ";
        cin >> file_name;

        proteinInputs = read_fasta_sequences(file_name);
        proteinIdInputs = read_fasta_ids(file_name);
        if (!proteinInputs.empty() && !proteinIdInputs.empty()) {
            break;
        }
    } while (true);

    do {
        file_name.clear();
        cout << "Input file name for DNA (exclude .fasta/.fastq):" << endl << "--> ";
        cin >> file_name;

        dnaInputs = read_fasta_sequences(file_name);
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

            int* t_sc_mat = (int*)malloc(size);
            int* t_ins_mat = (int*)malloc(size);
            int* t_del_mat = (int*)malloc(size);

            int* sc_mat_r = (int*)malloc(size);
            int* ins_mat_r = (int*)malloc(size);
            int* del_mat_r = (int*)malloc(size);

            int* t_sc_mat_r = (int*)malloc(size);
            int* t_ins_mat_r = (int*)malloc(size);
            int* t_del_mat_r = (int*)malloc(size);

            char* c_DNA_sequence = new char[N_size];
            char* c_protein_sequence = new char[M_size];
            char* c_DNA_sequence_r = new char[N_size];

            memcpy(c_DNA_sequence, DNA_sequence.c_str(), N_size);
            memcpy(c_protein_sequence, protein_sequence.c_str(), M_size);
            memcpy(c_DNA_sequence_r, DNA_sequence_r.c_str(), N_size);

            scoring_init(c_DNA_sequence, c_protein_sequence, sc_mat, ins_mat, del_mat, t_sc_mat,
                         t_ins_mat, t_del_mat, N, M);
            if (frame == 6)
                scoring_init(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, ins_mat_r, del_mat_r,
                             t_sc_mat_r, t_ins_mat_r, t_del_mat_r, N, M);

            if (mode == 0) {
                start = steady_clock::now();
                scoring_cpu(c_DNA_sequence, c_protein_sequence, sc_mat, ins_mat, del_mat, t_sc_mat,
                            t_ins_mat, t_del_mat, N, M);
                end = steady_clock::now();
                double elapsed1 = duration<double, milli>(end - start).count();

                traceV2_1d_check(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                                 index_prot, index);
                myArray[index_prot][3] = to_string(elapsed1);
                if (frame == 3)
                    top5(index[0], index_prot, index[1], index[2], top_scores, top_i, top_j,
                         top_indexes);
                // write_to_excel(index_dna, index_prot);
                total_r += elapsed1;
                cout << "Run DNA: " << index_dna << " Prot: " << index_prot << endl
                     << "Time in ms: " << elapsed1 << endl;
                cout << "Total runtime: " << total_r << endl;

                if (frame == 6) {
                    start = steady_clock::now();
                    scoring_cpu(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, ins_mat_r,
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
            } else if (mode == 1) {

                char* d_codon_aa;
                char* d_protein_sequence;
                char* d_codon_aa_r;

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

                checkCudaErrors(cudaMalloc(&d_codon_aa, N_size));
                checkCudaErrors(cudaMalloc(&d_protein_sequence, M_size));
                checkCudaErrors(cudaMalloc(&d_codon_aa_r, N_size));

                checkCudaErrors(cudaMemcpy(d_codon_aa,
                                           translate_each_position(c_DNA_sequence, N).data(),
                                           N_size, cudaMemcpyHostToDevice));
                checkCudaErrors(cudaMemcpy(d_protein_sequence, c_protein_sequence, M_size,
                                           cudaMemcpyHostToDevice));
                checkCudaErrors(cudaMemcpy(d_codon_aa_r,
                                           translate_each_position(c_DNA_sequence_r, N).data(),
                                           N_size, cudaMemcpyHostToDevice));

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
                        scoring_cuda<<<currGridDim, blockDimMain, 0, stream1>>>(
                            d_codon_aa, d_protein_sequence, d_sc_mat, d_ins_mat, d_del_mat,
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
                        scoring_cuda<<<currGridDim, blockDimMain, 0, stream1>>>(
                            d_codon_aa, d_protein_sequence, d_sc_mat, d_ins_mat, d_del_mat,
                            d_t_sc_mat, d_t_ins_mat, d_t_del_mat, (unsigned int)N, (unsigned int)M,
                            diag, submatrixSide, first_sub_y_curr_diag);
                        scoring_cuda<<<currGridDim, blockDimMain, 0, stream2>>>(
                            d_codon_aa_r, d_protein_sequence, d_sc_mat_r, d_ins_mat_r, d_del_mat_r,
                            d_t_sc_mat_r, d_t_ins_mat_r, d_t_del_mat_r, (unsigned int)N,
                            (unsigned int)M, diag, submatrixSide, first_sub_y_curr_diag);
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

                checkCudaErrors(cudaFree(d_codon_aa));
                checkCudaErrors(cudaFree(d_protein_sequence));
                checkCudaErrors(cudaFree(d_codon_aa_r));

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

            free(t_sc_mat);
            free(t_ins_mat);
            free(t_del_mat);

            free(sc_mat_r);
            free(ins_mat_r);
            free(del_mat_r);

            free(t_sc_mat_r);
            free(t_ins_mat_r);
            free(t_del_mat_r);

            delete[] c_protein_sequence;
            delete[] c_DNA_sequence;
            delete[] c_DNA_sequence_r;
        }

        for (int i = 0; i < 5; i++) {
            if (top == 0 && top_scores[i] != top_scores[0])
                continue;

            DNA_sequence = dnaInputs[index_dna];
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

            int* t_sc_mat = (int*)malloc(size);
            int* t_ins_mat = (int*)malloc(size);
            int* t_del_mat = (int*)malloc(size);

            int* sc_mat_r = (int*)malloc(size);
            int* ins_mat_r = (int*)malloc(size);
            int* del_mat_r = (int*)malloc(size);

            int* t_sc_mat_r = (int*)malloc(size);
            int* t_ins_mat_r = (int*)malloc(size);
            int* t_del_mat_r = (int*)malloc(size);

            char* c_DNA_sequence = new char[N_size];
            char* c_protein_sequence = new char[M_size];
            char* c_DNA_sequence_r = new char[N_size];

            memcpy(c_DNA_sequence, DNA_sequence.c_str(), N_size);
            memcpy(c_protein_sequence, protein_sequence.c_str(), M_size);
            memcpy(c_DNA_sequence_r, DNA_sequence_r.c_str(), N_size);

            scoring_init(c_DNA_sequence, c_protein_sequence, sc_mat, ins_mat, del_mat, t_sc_mat,
                         t_ins_mat, t_del_mat, N, M);
            if (frame == 6)
                scoring_init(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, ins_mat_r, del_mat_r,
                             t_sc_mat_r, t_ins_mat_r, t_del_mat_r, N, M);

            // use cpu here, the cuda code is so long and we only have to iterate another 5 times
            // its honestly worth it saving LOC and keeping code simple to run a couple ms slower
            cout << proteinIdInputs[top_indexes[i]] << endl
                 << "Sequence: " << top_indexes[i] << endl;
            start = steady_clock::now();
            scoring_cpu(c_DNA_sequence, c_protein_sequence, sc_mat, ins_mat, del_mat, t_sc_mat,
                        t_ins_mat, t_del_mat, N, M);
            end = steady_clock::now();
            double elapsed1 = duration<double, milli>(end - start).count();

            traceV2_1d_check(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                             top_indexes[i], index);

            if (frame == 3) {
                traceV2_1d(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                           top_indexes[i], index);
                cout << endl << "Score: " << top_scores[i] << endl;
                cout << "Start to End match in Protein: " << index[3] << "-" << index[2] << endl
                     << endl;
                cout << "Time in ms: " << elapsed1 << endl << endl;
            }

            if (frame == 6) {
                start = steady_clock::now();
                scoring_cpu(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, ins_mat_r, del_mat_r,
                            t_sc_mat_r, t_ins_mat_r, t_del_mat_r, N, M);
                end = steady_clock::now();
                double elapsed = duration<double, milli>(end - start).count();

                traceV2_1d_check(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, t_sc_mat_r, N, M,
                                 top_indexes[i], index_r);

                if (index[0] >= index_r[0]) {
                    traceV2_1d(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M,
                               top_indexes[i], index);
                    cout << endl << "Score: " << top_scores[i] << endl;
                    cout << "Start to End match in Protein: " << index[3] << "-" << index[2] << endl
                         << endl;
                    cout << "Time in ms: " << elapsed1 << endl << endl;
                } else {
                    cout << "Reverse: " << endl;
                    traceV2_1d(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, t_sc_mat_r, N, M,
                               top_indexes[i], index_r);
                    cout << endl << "Score: " << top_scores[i] << endl;
                    cout << "Start to End match in Protein: " << index_r[3] << "-" << index_r[2]
                         << endl
                         << endl;
                    cout << "Time in ms: " << elapsed << endl << endl;
                }
            }

            free(sc_mat);
            free(ins_mat);
            free(del_mat);

            free(t_sc_mat);
            free(t_ins_mat);
            free(t_del_mat);

            free(sc_mat_r);
            free(ins_mat_r);
            free(del_mat_r);

            free(t_sc_mat_r);
            free(t_ins_mat_r);
            free(t_del_mat_r);

            delete[] c_protein_sequence;
            delete[] c_DNA_sequence;
            delete[] c_DNA_sequence_r;
        }

        delete[] top_scores;
        delete[] top_i;
        delete[] top_j;
        delete[] top_indexes;
        delete[] index;
        delete[] index_r;
    }

    return 0;
}