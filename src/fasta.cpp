#include "fasta.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef CUDA_ALIGNMENT_DATA_DIR
#define CUDA_ALIGNMENT_DATA_DIR "./data"
#endif

using std::string;

std::vector<string> read_fasta_sequences(const string& filename) {
    std::ifstream file(string(CUDA_ALIGNMENT_DATA_DIR "/") + filename + ".fasta");
    std::vector<string> sequences;
    string line, sequence;

    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return sequences;
    }

    while (getline(file, line)) {
        if (line.empty())
            continue;
        if (line[0] == '>') {
            if (!sequence.empty()) {
                sequences.push_back(sequence);
                sequence.clear();
            }
        } else {
            sequence += line;
        }
    }

    if (!sequence.empty()) {
        sequences.push_back(sequence);
    }

    file.close();

    return sequences;
}

std::vector<string> read_fasta_ids(const string& filename) {
    std::ifstream file(string(CUDA_ALIGNMENT_DATA_DIR "/") + filename + ".fasta");
    std::vector<string> sequenceIDs;
    string line;

    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return sequenceIDs;
    }

    while (getline(file, line)) {
        if (!line.empty() && line[0] == '>') {
            sequenceIDs.push_back(line.substr(1)); // Remove '>' and store only the ID
        }
    }

    return sequenceIDs;
}
