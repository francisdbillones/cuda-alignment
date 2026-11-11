#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <algorithm>
#include <array>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <memory>

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

struct Hit {
    int score = 0;
    int i = 0, j = 0;
    bool reverse = false;
};

struct Ranked {
    int score = 0;
    int prot = 0;
};

std::array<Ranked, 5> top5(const vector<Hit>& hits) {
    std::array<Ranked, 5> top{};
    for (int p = 0; p < (int)hits.size(); p++) {
        for (int x = 0; x < 5; x++) {
            if (hits[p].score >= top[x].score) {
                for (int y = 4; y > x; y--)
                    top[y] = top[y - 1];
                top[x] = {hits[p].score, p};
                break;
            }
        }
    }
    return top;
}

struct Alignment {
    string dna, prot, frameshift;
};

void routine(int trace, int& i, int& j, char* str, char* ref_seq, Alignment& aln) {
    const int di = trace == TRACE_DEL ? 3 : trace == TRACE_INS ? 0 : trace;
    const int dj = trace != TRACE_DEL;
    aln.dna += di ? translate_codon(str + i - 1) : '-';
    aln.prot += dj ? ref_seq[j - 1] : '-';
    aln.frameshift += trace == TRACE_SHIFT_2 || trace == TRACE_SHIFT_4 ? '*' : ' ';
    i -= di;
    j -= dj;
}

void traceV2_1d_check(char* input_seq, char* ref_seq, int* sc_mat, int* t_sc_mat, size_t N,
                      size_t M, int* indexes, Alignment* aln = nullptr) {
    N = (int)N;
    M = (int)M;
    int i_max = 0, j_max = 0, max_score = 0, curr_score, i = 0, j = 0;
    Alignment out;

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

    while (sc_mat[i * M + j] != 0) {
        indexes[3] = j;
        routine(t_sc_mat[i * M + j], i, j, input_seq, ref_seq, out);
    }

    reverse(out.dna.begin(), out.dna.end());
    reverse(out.prot.begin(), out.prot.end());
    reverse(out.frameshift.begin(), out.frameshift.end());
    if (aln)
        *aln = out;
}

Hit best_cell(const int* sc_mat, size_t N, size_t M) {
    Hit h;
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            if (sc_mat[i * M + j] > h.score) {
                h.score = sc_mat[i * M + j];
                h.i = (int)i;
                h.j = (int)j;
            }
        }
    }
    return h;
}

Hit better_strand(Hit fwd, Hit rev) {
    if (fwd.score >= rev.score)
        return fwd;
    rev.reverse = true;
    return rev;
}

struct HostMats {
    std::unique_ptr<int[]> sc, ins, del, t_sc, t_ins, t_del;

    explicit HostMats(size_t n)
        : sc(new int[n]), ins(new int[n]), del(new int[n]), t_sc(new int[n]), t_ins(new int[n]),
          t_del(new int[n]) {}

    void init(char* dna, char* prot, size_t N, size_t M) {
        scoring_init(dna, prot, sc.get(), ins.get(), del.get(), t_sc.get(), t_ins.get(),
                     t_del.get(), N, M);
    }
};

struct DeviceMats {
    int *sc, *ins, *del, *t_sc, *t_ins, *t_del;

    DeviceMats(const HostMats& h, size_t n) {
        int** dst[] = {&sc, &ins, &del, &t_sc, &t_ins, &t_del};
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

vector<Hit> score_all_cpu(const string& dna, const vector<string>& prots, size_t num_prots,
                          int frame, double& elapsed_ms) {
    string fwd = dna, rev = reverse_complement(dna);
    const size_t N = dna.length();
    vector<Hit> hits(num_prots);

    for (size_t p = 0; p < num_prots; p++) {
        string prot = prots[p];
        const size_t M = prot.length() + 1;

        auto score_strand = [&](string& strand) {
            HostMats m(N * M);
            m.init(&strand[0], &prot[0], N, M);
            const auto start = steady_clock::now();
            scoring_cpu(strand.data(), prot.data(), m.sc.get(), m.ins.get(), m.del.get(),
                        m.t_sc.get(), m.t_ins.get(), m.t_del.get(), N, M);
            elapsed_ms += duration<double, milli>(steady_clock::now() - start).count();
            return best_cell(m.sc.get(), N, M);
        };

        hits[p] = score_strand(fwd);
        if (frame == 6)
            hits[p] = better_strand(hits[p], score_strand(rev));
    }
    return hits;
}

vector<Hit> score_all_cuda(const string& dna, const vector<string>& prots, size_t num_prots,
                           int frame, double& elapsed_ms) {
    string fwd = dna, rev = reverse_complement(dna);
    const size_t N = dna.length();
    vector<Hit> hits(num_prots);

    char *d_codon_aa, *d_codon_aa_r;
    checkCudaErrors(cudaMalloc(&d_codon_aa, N));
    checkCudaErrors(cudaMalloc(&d_codon_aa_r, N));
    checkCudaErrors(cudaMemcpy(d_codon_aa, translate_each_position(fwd.data(), N).data(), N,
                               cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_codon_aa_r, translate_each_position(rev.data(), N).data(), N,
                               cudaMemcpyHostToDevice));

    cudaStream_t stream1, stream2;
    checkCudaErrors(cudaStreamCreate(&stream1));
    checkCudaErrors(cudaStreamCreate(&stream2));
    GpuTimer timer;

    for (size_t p = 0; p < num_prots; p++) {
        string prot = prots[p];
        const size_t M = prot.length() + 1;
        const size_t n = N * M;

        HostMats h(n);
        h.init(&fwd[0], &prot[0], N, M);
        std::unique_ptr<HostMats> h_r;
        if (frame == 6) {
            h_r = std::make_unique<HostMats>(n);
            h_r->init(&rev[0], &prot[0], N, M);
        }

        char* d_protein_sequence;
        checkCudaErrors(cudaMalloc(&d_protein_sequence, M));
        checkCudaErrors(cudaMemcpy(d_protein_sequence, prot.data(), M, cudaMemcpyHostToDevice));
        DeviceMats d(h, n);
        std::unique_ptr<DeviceMats> d_r;
        if (frame == 6)
            d_r = std::make_unique<DeviceMats>(*h_r, n);

        dim3 blockDimMain(32, 32);
        unsigned int submatrixSide = blockDimMain.y;
        unsigned int numSubmatrixRows = ((unsigned int)N + submatrixSide - 1) / submatrixSide;
        unsigned int numSubmatrixCols = ((unsigned int)M + submatrixSide - 1) / submatrixSide;

        timer.Start();
        for (unsigned int diag = 0; diag < numSubmatrixRows + numSubmatrixCols - 1; ++diag) {
            unsigned int first_sub_y_curr_diag =
                std::max(0, (int)diag - (int)(numSubmatrixCols - 1));
            unsigned int last_sub_y_curr_diag = std::min((int)diag, (int)numSubmatrixRows - 1);
            dim3 currGridDim(last_sub_y_curr_diag - first_sub_y_curr_diag + 1, 1);

            scoring_cuda<<<currGridDim, blockDimMain, 0, stream1>>>(
                d_codon_aa, d_protein_sequence, d.sc, d.ins, d.del, d.t_sc, d.t_ins, d.t_del,
                (unsigned int)N, (unsigned int)M, diag, submatrixSide, first_sub_y_curr_diag);
            if (d_r)
                scoring_cuda<<<currGridDim, blockDimMain, 0, stream2>>>(
                    d_codon_aa_r, d_protein_sequence, d_r->sc, d_r->ins, d_r->del, d_r->t_sc,
                    d_r->t_ins, d_r->t_del, (unsigned int)N, (unsigned int)M, diag, submatrixSide,
                    first_sub_y_curr_diag);
        }
        checkCudaErrors(cudaStreamSynchronize(stream1));
        checkCudaErrors(cudaStreamSynchronize(stream2));
        timer.Stop();
        elapsed_ms += timer.Elapsed();

        checkCudaErrors(cudaMemcpy(h.sc.get(), d.sc, n * sizeof(int), cudaMemcpyDeviceToHost));
        hits[p] = best_cell(h.sc.get(), N, M);
        if (d_r) {
            checkCudaErrors(
                cudaMemcpy(h_r->sc.get(), d_r->sc, n * sizeof(int), cudaMemcpyDeviceToHost));
            hits[p] = better_strand(hits[p], best_cell(h_r->sc.get(), N, M));
        }
        checkCudaErrors(cudaFree(d_protein_sequence));
    }

    checkCudaErrors(cudaStreamDestroy(stream1));
    checkCudaErrors(cudaStreamDestroy(stream2));
    checkCudaErrors(cudaFree(d_codon_aa));
    checkCudaErrors(cudaFree(d_codon_aa_r));
    return hits;
}

int main() {
    int mode, frame, top, device;
    vector<string> dnaInputs, proteinInputs, proteinIdInputs;
    string protein_sequence, DNA_sequence, DNA_sequence_r, file_name;

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
    const size_t num_prots = min<size_t>(proteinInputs.size(), 4096);

    // HITS_DUMP=path dumps
    // "read protein score i j reverse"
    // for everything
    // need this to make sure our kernel still generates the same thing
    std::ofstream hits_dump;
    if (const char* path = std::getenv("HITS_DUMP"))
        hits_dump.open(path);

    for (int index_dna = 0; index_dna < dnaInputs.size(); index_dna++) {
        const vector<Hit> hits =
            mode == 0
                ? score_all_cpu(dnaInputs[index_dna], proteinInputs, num_prots, frame, total_r)
                : score_all_cuda(dnaInputs[index_dna], proteinInputs, num_prots, frame, total_r);
        cout << "Total runtime: " << total_r << endl;

        if (hits_dump.is_open())
            for (size_t p = 0; p < hits.size(); p++)
                hits_dump << index_dna << ' ' << p << ' ' << hits[p].score << ' ' << hits[p].i
                          << ' ' << hits[p].j << ' ' << hits[p].reverse << '\n';

        const std::array<Ranked, 5> best = top5(hits);
        int index[4] = {}, index_r[4] = {};
        steady_clock::time_point start, end;

        for (int i = 0; i < 5; i++) {
            if (top == 0 && best[i].score != best[0].score)
                continue;

            DNA_sequence = dnaInputs[index_dna];
            protein_sequence = proteinInputs[best[i].prot];
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

            char* c_DNA_sequence = new char[N_size + 1];
            char* c_protein_sequence = new char[M_size];
            char* c_DNA_sequence_r = new char[N_size + 1];

            memcpy(c_DNA_sequence, DNA_sequence.c_str(), N_size + 1);
            memcpy(c_protein_sequence, protein_sequence.c_str(), M_size);
            memcpy(c_DNA_sequence_r, DNA_sequence_r.c_str(), N_size + 1);

            scoring_init(c_DNA_sequence, c_protein_sequence, sc_mat, ins_mat, del_mat, t_sc_mat,
                         t_ins_mat, t_del_mat, N, M);
            if (frame == 6)
                scoring_init(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, ins_mat_r, del_mat_r,
                             t_sc_mat_r, t_ins_mat_r, t_del_mat_r, N, M);

            // use cpu here, the cuda code is so long and we only have to iterate another 5 times
            // its honestly worth it saving LOC and keeping code simple to run a couple ms slower
            cout << proteinIdInputs[best[i].prot] << endl << "Sequence: " << best[i].prot << endl;
            start = steady_clock::now();
            scoring_cpu(c_DNA_sequence, c_protein_sequence, sc_mat, ins_mat, del_mat, t_sc_mat,
                        t_ins_mat, t_del_mat, N, M);
            end = steady_clock::now();
            double elapsed1 = duration<double, milli>(end - start).count();

            Alignment aln, aln_r;
            traceV2_1d_check(c_DNA_sequence, c_protein_sequence, sc_mat, t_sc_mat, N, M, index,
                             &aln);

            bool use_reverse = false;
            double elapsed = elapsed1;
            if (frame == 6) {
                start = steady_clock::now();
                scoring_cpu(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, ins_mat_r, del_mat_r,
                            t_sc_mat_r, t_ins_mat_r, t_del_mat_r, N, M);
                end = steady_clock::now();
                double elapsed_r = duration<double, milli>(end - start).count();

                traceV2_1d_check(c_DNA_sequence_r, c_protein_sequence, sc_mat_r, t_sc_mat_r, N, M,
                                 index_r, &aln_r);
                if (index[0] < index_r[0]) {
                    use_reverse = true;
                    elapsed = elapsed_r;
                }
            }

            const char* dna = use_reverse ? c_DNA_sequence_r : c_DNA_sequence;
            Alignment& shown = use_reverse ? aln_r : aln;
            const int* shown_index = use_reverse ? index_r : index;
            string f1, f2, f3;
            size_t k;

            if (use_reverse)
                cout << "Reverse: " << endl;
            cout << endl;

            three_frame(dna, &f1, &f2, &f3);

            cout << "DNA Input:\t";
            for (k = 0; k < strlen(dna); k++) {
                cout << dna[k];
            }
            cout << endl;
            cout << "\t\t ";
            for (k = 0; k < f1.length(); k++) {
                cout << f1[k] << "  ";
            }
            cout << endl;
            cout << "\t\t  ";
            for (k = 0; k < f2.length(); k++) {
                cout << f2[k] << "  ";
            }
            cout << endl;
            cout << "\t\t   ";
            for (k = 0; k < f3.length(); k++) {
                cout << f3[k] << "  ";
            }
            cout << endl << endl;

            cout << "Frame Match: \t ";
            for (k = 0; k < shown.dna.length(); k++) {
                cout << shown.dna[k] << "  ";
            }
            cout << endl;
            cout << "Frameshift: \t ";
            for (k = 0; k < shown.frameshift.length(); k++) {
                cout << shown.frameshift[k] << "  ";
            }
            cout << endl;
            cout << "Output Prot: \t ";
            for (k = 0; k < shown.prot.length(); k++) {
                cout << shown.prot[k] << "  ";
            }
            cout << endl;

            cout << endl << "Score: " << best[i].score << endl;
            cout << "Start to End match in Protein: " << shown_index[3] << "-" << shown_index[2]
                 << endl
                 << endl;
            cout << "Time in ms: " << elapsed << endl << endl;

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
    }

    return 0;
}