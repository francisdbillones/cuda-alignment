#ifndef BLOSUM62_H
#define BLOSUM62_H

extern int blosum62mat[24][24];

extern int prot_to_idx[128];

void read_blosum62();

void initialize_protein_to_int_map();

#endif // BLOSUM62_H
