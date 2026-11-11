#ifndef CLI_H
#define CLI_H

#include <string>
#include <vector>

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
               const Alignment& aln, int score, const int* index, double ms);

#endif // CLI_H
