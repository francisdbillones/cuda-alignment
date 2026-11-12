#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include <algorithm>
#include <array>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <memory>

#include "blosum62.h"
#include "cli.h"
#include "codon.cuh"
#include "cuda_check.cuh"
#include "fasta.h"
#include "gputimer.h"
#include "mats.cuh"
#include "scoring.cuh"

using std::string;
using std::vector;

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

void routine(int trace, int& i, int& j, char* str, char* ref_seq, Alignment& aln) {
    const int di = trace == TRACE_DEL ? 3 : trace == TRACE_INS ? 0 : trace;
    const int dj = trace != TRACE_DEL;
    aln.dna += di ? translate_codon(str + i - 1) : '-';
    aln.prot += dj ? ref_seq[j - 1] : '-';
    aln.frameshift += trace == TRACE_SHIFT_2 || trace == TRACE_SHIFT_4 ? '*' : ' ';
    i -= di;
    j -= dj;
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

struct Traceback {
    int start_j = 0, end_j = 0;
    Alignment aln;
};

Traceback traceV2_1d_check(char* input_seq, char* ref_seq, const HostMats& m, size_t N, size_t M) {
    const int* sc_mat = m.sc.get();
    const int* t_sc_mat = m.t_sc.get();
    const Hit best = best_cell(sc_mat, N, M);
    int i = best.i, j = best.j;
    Traceback t;
    t.end_j = j;

    while (sc_mat[i * M + j] != 0) {
        t.start_j = j;
        routine(t_sc_mat[i * M + j], i, j, input_seq, ref_seq, t.aln);
    }

    reverse(t.aln.dna.begin(), t.aln.dna.end());
    reverse(t.aln.prot.begin(), t.aln.prot.end());
    reverse(t.aln.frameshift.begin(), t.aln.frameshift.end());
    return t;
}

Hit better_strand(Hit fwd, Hit rev) {
    if (fwd.score >= rev.score)
        return fwd;
    rev.reverse = true;
    return rev;
}

double run_cpu(string& dna, string& prot, HostMats& m) {
    const size_t N = dna.length(), M = prot.length() + 1;
    m.init(&dna[0], &prot[0], N, M);
    const auto start = std::chrono::steady_clock::now();
    scoring_cpu(dna.data(), prot.data(), m.sc.get(), m.ins.get(), m.del.get(), m.t_sc.get(),
                m.t_ins.get(), m.t_del.get(), N, M);
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
        .count();
}

Traceback trace_cpu(string& dna, string& prot, double& ms) {
    const size_t N = dna.length(), M = prot.length() + 1;
    HostMats m(N * M);
    ms = run_cpu(dna, prot, m);
    return traceV2_1d_check(&dna[0], &prot[0], m, N, M);
}

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
            elapsed_ms += run_cpu(strand, prot, m);
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
    int device;
    checkCudaErrors(cudaGetDevice(&device));
    checkCudaErrors(cudaSetDevice(device));

    read_blosum62();
    initialize_protein_to_int_map();
    upload_blosum62_tables();

    const Options opt = prompt_options();
    double total_r = 0.0;
    const size_t num_prots = std::min<size_t>(opt.prots.size(), 4096);

    // HITS_DUMP=path dumps
    // "read protein score i j reverse"
    // for everything
    // need this to make sure our kernel still generates the same thing
    std::ofstream hits_dump;
    if (const char* path = std::getenv("HITS_DUMP"))
        hits_dump.open(path);

    for (int index_dna = 0; index_dna < opt.dnas.size(); index_dna++) {
        const vector<Hit> hits =
            opt.mode == 0
                ? score_all_cpu(opt.dnas[index_dna], opt.prots, num_prots, opt.frame, total_r)
                : score_all_cuda(opt.dnas[index_dna], opt.prots, num_prots, opt.frame, total_r);
        print_runtime(total_r);

        if (hits_dump.is_open())
            for (size_t p = 0; p < hits.size(); p++)
                hits_dump << index_dna << ' ' << p << ' ' << hits[p].score << ' ' << hits[p].i
                          << ' ' << hits[p].j << ' ' << hits[p].reverse << '\n';

        const std::array<Ranked, 5> best = top5(hits);
        string fwd = opt.dnas[index_dna], rev = reverse_complement(fwd);

        for (int i = 0; i < 5; i++) {
            if (opt.top == 0 && best[i].score != best[0].score)
                continue;

            string prot = opt.prots[best[i].prot];
            const bool use_reverse = hits[best[i].prot].reverse;
            string& dna = use_reverse ? rev : fwd;

            // use cpu here, the cuda code is so long and we only have to iterate another 5 times
            // its honestly worth it saving LOC and keeping code simple to run a couple ms slower
            double elapsed;
            const Traceback t = trace_cpu(dna, prot, elapsed);
            print_hit(opt.prot_ids[best[i].prot], best[i].prot, dna, use_reverse, t.aln,
                      best[i].score, t.start_j, t.end_j, elapsed);
        }
    }

    return 0;
}