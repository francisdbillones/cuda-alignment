#ifndef CODON_CUH
#define CODON_CUH

#include <algorithm>
#include <string>

__host__ __device__ inline int nucleotide_index(const char base) {
    switch (base) {
    case 'T':
        return 0;
    case 'C':
        return 1;
    case 'A':
        return 2;
    case 'G':
        return 3;
    default:
        return -1;
    }
}

__host__ __device__ inline char translate_codon(const char* codon) {
    const char* table = "FFLL"
                        "SSSS"
                        "YYXX"
                        "CCXW"
                        "LLLL"
                        "PPPP"
                        "HHQQ"
                        "RRRR"
                        "IIIM"
                        "TTTT"
                        "NNKK"
                        "SSRR"
                        "VVVV"
                        "AAAA"
                        "DDEE"
                        "GGGG";

    const int first = nucleotide_index(codon[0]);
    const int second = nucleotide_index(codon[1]);
    const int third = nucleotide_index(codon[2]);
    if (first < 0 || second < 0)
        return ' ';

    const char* box = table + 16 * first + 4 * second;
    if (third >= 0)
        return box[third];

    // if all characters match then third base irrelevant
    if (box[0] == box[1] && box[1] == box[2] && box[2] == box[3])
        return box[0];
    else
        return ' ';
}

__host__ __device__ inline char complement(const char base) {
    switch (base) {
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'C':
        return 'G';
    case 'G':
        return 'C';
    default:
        return '\0';
    }
}

inline std::string reverse_complement(std::string dna) {
    std::reverse(dna.begin(), dna.end());
    std::transform(dna.begin(), dna.end(), dna.begin(), complement);
    return dna;
}

#endif // CODON_CUH
