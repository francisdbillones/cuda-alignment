#include "cli.h"

#include <iostream>

#include "codon.cuh"
#include "fasta.h"

using std::cout;
using std::endl;
using std::string;

static string prompt_file(const char* msg) {
    string file_name;
    cout << msg << endl << "--> ";
    std::cin >> file_name;
    return file_name;
}

static int prompt_choice(const char* msg, int a, int b) {
    while (true) {
        int choice;
        cout << msg << endl << "--> ";
        std::cin >> choice;
        if (choice == a || choice == b)
            return choice;
        cout << "Invalid input please try again." << endl;
    }
}

Options prompt_options() {
    Options o;
    do {
        const string f = prompt_file("Input file name for protein(exclude .fasta/.fastq):");
        o.prots = read_fasta_sequences(f);
        o.prot_ids = read_fasta_ids(f);
    } while (o.prots.empty() || o.prot_ids.empty());
    do {
        o.dnas =
            read_fasta_sequences(prompt_file("Input file name for DNA (exclude .fasta/.fastq):"));
    } while (o.dnas.empty());
    o.mode = prompt_choice("Choose mode [0 - Sequential] [1 - CUDA]:", 0, 1);
    o.frame = prompt_choice("Choose frame count [3 - 3 Frame] [6 - 6 Frame]:", 3, 6);
    o.top = prompt_choice("Show Score [0 - Top 1] [1 - Top 5]:", 0, 1);
    return o;
}

void print_runtime(double ms) {
    cout << "Total runtime: " << ms << endl;
}

static void print_spaced(const char* label, const string& s) {
    cout << label;
    for (char c : s)
        cout << c << "  ";
    cout << endl;
}

void print_hit(const string& id, int prot, const string& dna, bool reverse, const Alignment& aln,
               int score, int start_j, int end_j, double ms) {
    cout << id << endl << "Sequence: " << prot << endl;
    if (reverse)
        cout << "Reverse: " << endl;
    cout << endl;

    string f1, f2, f3;
    three_frame(dna.c_str(), &f1, &f2, &f3);
    cout << "DNA Input:\t" << dna << endl;
    print_spaced("\t\t ", f1);
    print_spaced("\t\t  ", f2);
    print_spaced("\t\t   ", f3);
    cout << endl;
    print_spaced("Frame Match: \t ", aln.dna);
    print_spaced("Frameshift: \t ", aln.frameshift);
    print_spaced("Output Prot: \t ", aln.prot);

    cout << endl << "Score: " << score << endl;
    cout << "Start to End match in Protein: " << start_j << "-" << end_j << endl << endl;
    cout << "Time in ms: " << ms << endl << endl;
}
