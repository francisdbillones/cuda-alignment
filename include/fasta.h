#ifndef FASTA_H
#define FASTA_H

#include <string>
#include <vector>

// Both read CUDA_ALIGNMENT_DATA_DIR/<filename>.fasta. On open failure they print
// an error and return an empty vector.
std::vector<std::string> readFastaSequences(const std::string& filename);
std::vector<std::string> readFastaIDs(const std::string& filename);

#endif // FASTA_H
