#ifndef CLI_H
#define CLI_H

#include <string>
#include <vector>

constexpr int TOP_N = 5;

struct Alignment {
    std::string dna, prot, frameshift;
};

struct Options {
    std::vector<std::string> prots, prot_ids, dnas;
    int mode, frame, top;
};

Options prompt_options();

void print_runtime(double ms);

void print_hit(const std::string& id, int prot, const std::string& dna, bool reverse,
               const Alignment& aln, int score, int start_j, int end_j, double ms);

#endif // CLI_H
