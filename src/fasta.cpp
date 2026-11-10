#include "fasta.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef CUDA_ALIGNMENT_DATA_DIR
#define CUDA_ALIGNMENT_DATA_DIR "./data"
#endif

using std::cerr;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;

vector<string> readFastaSequences(const string& filename) {
    ifstream file(string(CUDA_ALIGNMENT_DATA_DIR "/") + filename + ".fasta");
    vector<string> sequences;
    string line, sequence;

    if (!file) {
        cerr << "Error: Unable to open file " << filename << endl;
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

vector<string> readFastaIDs(const string& filename) {
    ifstream file(string(CUDA_ALIGNMENT_DATA_DIR "/") + filename + ".fasta");
    vector<string> sequenceIDs;
    string line;

    if (!file) {
        cerr << "Error: Unable to open file " << filename << endl;
        return sequenceIDs;
    }

    while (getline(file, line)) {
        if (!line.empty() && line[0] == '>') {
            sequenceIDs.push_back(line.substr(1)); // Remove '>' and store only the ID
        }
    }

    return sequenceIDs;
}
