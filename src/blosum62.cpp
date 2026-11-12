#include "blosum62.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef CUDA_ALIGNMENT_DATA_DIR
#define CUDA_ALIGNMENT_DATA_DIR "./data"
#endif

int blosum62mat[24][24];

int prot_to_idx[128];

void initialize_protein_to_int_map() {

    std::vector<int> host_map(128, 22);

    host_map['A'] = 0;
    host_map['R'] = 1;
    host_map['N'] = 2;
    host_map['D'] = 3;
    host_map['C'] = 4;
    host_map['Q'] = 5;
    host_map['E'] = 6;
    host_map['G'] = 7;
    host_map['H'] = 8;
    host_map['I'] = 9;
    host_map['L'] = 10;
    host_map['K'] = 11;
    host_map['M'] = 12;
    host_map['F'] = 13;
    host_map['P'] = 14;
    host_map['S'] = 15;
    host_map['T'] = 16;
    host_map['W'] = 17;
    host_map['Y'] = 18;
    host_map['V'] = 19;
    host_map['B'] = 20;
    host_map['Z'] = 21;
    host_map['X'] = 22;
    host_map['U'] = 23;

    for (int i = 0; i < 128; ++i)
        prot_to_idx[i] = host_map[i];
}

void read_blosum62() {
    std::ifstream file(CUDA_ALIGNMENT_DATA_DIR "/BLOSUM62.txt");
    if (!file.is_open()) {
        std::cerr << "Error opening the BLOSUM62 file!" << std::endl;
        exit(1);
    }
    std::string skip;
    std::getline(file, skip);
    std::getline(file, skip);

    char c = 0;

    for (auto& row : blosum62mat) {
        file >> c;
        for (int& entry : row) {
            file >> entry;
        }
    }
    file.close();
}
