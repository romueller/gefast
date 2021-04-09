/*
 * GeFaST
 *
 * Copyright (C) 2016 - 2021 Robert Mueller
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact: Robert Mueller <romueller@techfak.uni-bielefeld.de>
 * Faculty of Technology, Bielefeld University,
 * PO box 100131, DE-33501 Bielefeld, Germany
 */

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>

#include "../libs/stats/include/stats.hpp"

#include "../include/Utility.hpp"

namespace GeFaST {

    void print_information() {

        std::cout << "##### GeFaST (" + VERSION + ") #####" << std::endl;
        std::cout << "Copyright (C) 2016 - 2021 Robert Mueller" << std::endl;
        std::cout << "https://github.com/romueller/gefast" << std::endl << std::endl;

    }


    void print_help() {

        std::cout << "Usage: GeFaST <mode> <input> <config> [option]..." << std::endl << std::endl;

        std::cout << "The first three arguments are mandatory and fixed in their order:" << std::endl;
        std::cout << " <mode>     The abbreviation of the mode. " << std::endl;
        std::cout << " <input>    By default, GeFaST expects a comma-separated list of file paths as its second argument." << std::endl;
        std::cout << "            This behaviour can be changed with the --list_file option (see input / output options below)." << std::endl;
        std::cout << " <config>   The file path of the configuration file containing the (basic) configuration for this execution of GeFaST." << std::endl;
        std::cout << "These arguments are followed by an arbitrary list of options described below." << std::endl << std::endl;

        std::cout << "Modes:" << std::endl;
        std::cout << " lev     Cluster amplicons based on the number of edit operations in optimal pairwise alignments." << std::endl;
        std::cout << " as      Cluster amplicons based on the score of optimal pairwise alignments." << std::endl;
        std::cout << " qlev    Cluster amplicons based on the number of edit operations in optimal pairwise alignments considering the quality scores associated with the sequences." << std::endl;
        std::cout << " qas     Cluster amplicons based on the score of optimal pairwise alignments considering the quality scores associated with the sequences" << std::endl;
        std::cout << " cons    Cluster amplicons using a notion of consistency considering the quality and abundance of amplicons." << std::endl;
        std::cout << " derep   Group amplicons based on exact sequence equality" << std::endl << std::endl;

        std::cout << "Configuration file:" << std::endl;
        std::cout << "The configuration file is a (possibly empty) plain-text file containing key-value pairs." << std::endl;
        std::cout << "Each configuration parameter is written in its own line and has the form <key>=<value>." << std::endl;
        std::cout << "Empty lines and comment lines (starting with #) are allowed." << std::endl;
        std::cout << "When an option is specified in the configuration file and on the command line, the value provided on the command line is used." << std::endl << std::endl;

        std::cout << "Input / output options:" << std::endl;
        std::cout << " -a, --alphabet STRING                 discard sequences with other characters (default: ACGT)" << std::endl;
        std::cout << " -i, --output_internal STRING          output links underlying the cluster to file (default: not created)" << std::endl;
        std::cout << " -o, --output_otus STRING              output clusters to file (default: not created)" << std::endl;
        std::cout << " -s, --output_statistics STRING        output statistics to file (defaut: not created)" << std::endl;
        std::cout << " -w, --output_seeds STRING             output seeds to file (default: not created)" << std::endl;
        std::cout << " -u, --output_uclust STRING            create UCLUST-like output file (default: not created)" << std::endl;
        std::cout << " --sep_abundance STRING                change separator symbol between identifier and abundance (default: _)" << std::endl;
        std::cout << " --min_length INTEGER                  discard shorter sequences (default: deactivated)" << std::endl;
        std::cout << " --max_length INTEGER                  discard longer sequences (default: deactivated)" << std::endl;
        std::cout << " --min_abundance INTEGER               discard less abundant sequences (default: deactivated)" << std::endl;
        std::cout << " --max_abundance INTEGER               discard more abundant sequences (default: deactivated)" << std::endl;
        std::cout << " --mothur                              output clusters in mothur list file format" << std::endl;
        std::cout << " --quality_encoding STRING             change expected quality encoding (FASTQ inputs, default: sanger)" << std::endl;
        std::cout << " --list_file STRING                    consider <input> option as path to file containing list of actual input files (one path per line)" << std::endl << std::endl;

        std::cout << "Clustering / refinement options:" << std::endl;
        std::cout << " -t, --threshold NUMERIC               distance threshold in clustering phase (default: mode-dependent)" << std::endl;
        std::cout << " -r, --refinement_threshold NUMERIC    distance threshold in refinement phase (default: 0, i.e. no refinement)" << std::endl;
        std::cout << " -b, --boundary INTEGER                mass boundary distinguishing between light and heavy clusters during refinement (default: 3)" << std::endl;
        std::cout << " -n, --break_swarms BINARY             do not extend cluster when the new amplicon has a larger abundance than the current subseed (default: 1)" << std::endl;
        std::cout << " -m, --match_reward INTEGER            reward for nucleotide match during pairwise global alignment computation (default: 5)" << std::endl;
        std::cout << " -p, --mismatch_penalty INTEGER        penalty for nucleotide mismatch during pairwise global alignment computation (default: -4)" << std::endl;
        std::cout << " -g, --gap_opening_penalty INTEGER     penalty for opening a gap during pairwise global alignment computation (default: -12)" << std::endl;
        std::cout << " -e, --gap_extension_penalty INTEGER   penalty for extending a gap during pairwise global alignment computation (default: -4)" << std::endl << std::endl;

        std::cout << "Component options:" << std::endl;
        std::cout << " --preprocessor STRING                 use the specified component to perform the preprocessing (default: mode-dependent)" << std::endl;
        std::cout << " --clusterer STRING                    use the specified component to cluster the amplicons (default: mode-dependent)" << std::endl;
        std::cout << " --refiner STRING                      use the specified component to refine the clusters (default: mode-dependent)" << std::endl;
        std::cout << " --output_generator STRING             use the specified component to generate the requested outputs (default: mode-dependent)" << std::endl << std::endl;

        std::cout << "See the manual for more details and mode-specific options." << std::endl;

    }


    std::vector<std::string> read_list_file(const std::string& list_file) {

        std::vector<std::string> files;

        std::ifstream in_stream(list_file);
        if (!in_stream.good()) {

            std::cerr << "ERROR: List file '" << list_file << "' not opened correctly." << std::endl;
            return files;

        }

        std::string file;

        while (std::getline(in_stream, file).good()) {

            if (file.empty() || file[0] == '#') continue; // skip empty and comment lines (begin with '#')

            files.push_back(file);

        }

        return files;

    }


    bool check_sequence(const std::string& seq, const std::string& alphabet, numSeqs_t seq_abund, lenSeqs_t min_len,
            lenSeqs_t max_len, numSeqs_t min_abund, numSeqs_t max_abund) {

        bool valid = true;

        if (!alphabet.empty()) {
            valid = seq.find_first_not_of(alphabet) == std::string::npos;
        }

        valid = valid
                && min_len <= seq.length() && seq.length() <= max_len
                && min_abund <= seq_abund && seq_abund <= max_abund;

        return valid;

    }


    char convert[128] = { // upper-case for a-z, everything else remains the same
            0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15,
            16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
            32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
            48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
            64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
            80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,
            96,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
            80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90, 123, 124, 125, 126, 127
    };

    void upper_case(std::string& s) {

        //std::transform(s.begin(), s.end(), s.begin(), ::toupper); // safer, more general (but slower) version
        for (auto i = 0; i < s.size(); i++) {
            s[i] = convert[s[i]];
        }

    }


    unsigned long long gcd(unsigned long long a, unsigned long long b) {
        return (b == 0) ? a : gcd(b, a % b);
    }


    /* === Defline === */

    Defline::Defline(const std::string& i, const std::string& info, const numSeqs_t a) :
        id(i), information(info), abundance(a) {
        // nothing else to do
    }

    Defline Defline::parse_description_line(const std::string& defline, const std::string& sep) {

        auto pos = defline.find(' ');
        numSeqs_t abundance = 1;
        std::string first_part = defline.substr(1, pos - 1);
        std::string second_part = defline.substr(pos + 1);

        pos = first_part.find(sep);

        if (pos != std::string::npos) {
            abundance = std::stoul(first_part.substr(pos + sep.size()));
            first_part = first_part.substr(0, pos);
        }

        return Defline(first_part, second_part, abundance);

    }


    /* === StringIteratorPair === */

    //TODO? custom hash function (e.g. FNV) to avoid construction of temporary string
    size_t hashStringIteratorPair::operator()(const StringIteratorPair& p) const {
        return hash(std::string(p.first, p.second));
    }

    bool equalStringIteratorPair::operator()(const StringIteratorPair& lhs, const StringIteratorPair& rhs) const {
        return ((lhs.second - lhs.first) == (rhs.second - rhs.first)) && std::equal(lhs.first, lhs.second, rhs.first);
    }

    bool lessStringIteratorPair::operator()(const StringIteratorPair& a, const StringIteratorPair& b) const {
        return std::lexicographical_compare(a.first, a.second, b.first, b.second);
    }


    /* === StdPair === */

    size_t hashStdPair::operator()(const std::pair<lenSeqs_t, lenSeqs_t>& p) const {

        size_t seed = hash(p.first);
        return hash(p.second) + 0x9e3779b9 + (seed << 6ul) + (seed >> 2ul);

    }


    /* === Dada2Utility === */

    void Dada2Utility::nt2int(char* oseq, const char* iseq) {

        int len = strlen(iseq);

        for (int i = 0; i < len; i++, iseq++, oseq++) {
            switch (*iseq) {
                case 'A':
                    *oseq = 1;
                    break;
                case 'C':
                    *oseq = 2;
                    break;
                case 'G':
                    *oseq = 3;
                    break;
                case 'T':
                    *oseq = 4;
                    break;
                case 'N':
                    *oseq = 5;
                    break;
                case '-':
                    *oseq = '-';
                    break;
                default:
                    printf("invalid character in input:%c.\n", *iseq);
                    *oseq = '\0';
            }
        }
        *oseq = '\0';

    }

    void Dada2Utility::int2nt(char* oseq, const char* iseq) {

        int len = strlen(iseq);

        for (int i = 0; i < len; i++, iseq++, oseq++) {
            switch (*iseq) {
                case 1:
                    *oseq = 'A';
                    break;
                case 2:
                    *oseq = 'C';
                    break;
                case 3:
                    *oseq = 'G';
                    break;
                case 4:
                    *oseq = 'T';
                    break;
                case 5:
                    *oseq = 'N';
                    break;
                case '-':
                    *oseq = '-';
                    break;
                default:
                    break;
            }
        }
        *oseq = '\0';

    }


    void Dada2Utility::assign_kmer8(uint8_t* kvec8, const char* seq, int k) {  // assumes a clean sequence (just 1s, 2s, 3s, 4s)

        size_t len = strlen(seq);
        if (len <= 0) throw std::runtime_error("Unexpected sequence length.");
        if (k >= len || k < 3 || k > 8) throw std::runtime_error("Invalid k-mer size.");

        size_t klen = len - k + 1; // number of k-mers in this sequence
        size_t n_kmers = (1 << (2 * k));  // 4^k k-mers
        uint16_t* kvec = new uint16_t[n_kmers];
        for (size_t kmer = 0; kmer < n_kmers; kmer++) {
            kvec[kmer] = 0;
        }

        for (int i = 0; i < klen; i++) {

            size_t kmer = 0;
            for (int j = i; j < i + k; j++) {

                int nti = ((int)seq[j]) - 1; // change 1s, 2s, 3s, 4s, to 0/1/2/3
                if (nti != 0 && nti != 1 && nti != 2 && nti != 3) {
                    throw std::runtime_error("Unexpected nucleotide.");
                }
                kmer = 4 * kmer + nti;

            }

            // Make sure k-mer index is valid. This doesn't solve the N's/-'s
            // issue though, as the "length" of the string (# of k-mers) needs
            // to also reflect the reduction from the N's/-'s
            if (kmer == 999999) {
                ;
            } else if (kmer >= n_kmers) {
                throw std::runtime_error("Kmer index out of range.");
            } else { // valid k-mer
                kvec[kmer]++;
            }

        }

        // move to destination uint8_t vector w/ overflow checks
        for (size_t kmer = 0; kmer < n_kmers; kmer++) {
            if (kvec[kmer]<255) {
                kvec8[kmer] = (uint8_t)kvec[kmer];
            } else {
                kvec8[kmer] = 255;
            }
        }

        delete[] kvec;

    }

    void Dada2Utility::assign_kmer(uint16_t* kvec, const char* seq, int k) {  // assumes a clean sequence (just 1s, 2s, 3s, 4s)

        size_t len = strlen(seq);
        if (len <= 0) throw std::runtime_error("Unexpected sequence length.");
        if (k >= len || k < 3 || k > 8) throw std::runtime_error("Invalid k-mer size.");

        size_t klen = len - k + 1; // number of k-mers in this sequence
        size_t n_kmers = (1 << (2 * k));  // 4^k k-mers
        for (size_t kmer = 0; kmer < n_kmers; kmer++) {
            kvec[kmer] = 0;
        }

        for (int i = 0; i < klen; i++) {

            size_t kmer = 0;
            for (int j = i; j < i + k; j++) {

                int nti = ((int)seq[j]) - 1; // change 1s, 2s, 3s, 4s, to 0/1/2/3
                if (nti != 0 && nti != 1 && nti != 2 && nti != 3) {
                    throw std::runtime_error("Unexpected nucleotide.");
                }
                kmer = 4 * kmer + nti;

            }

            // Make sure k-mer index is valid. This doesn't solve the N's/-'s
            // issue though, as the "length" of the string (# of k-mers) needs
            // to also reflect the reduction from the N's/-'s
            if (kmer == 999999) {
                ;
            } else if (kmer >= n_kmers) {
                throw std::runtime_error("Kmer index out of range.");
            } else { // valid k-mer
                kvec[kmer]++;
            }

        }

    }

    void Dada2Utility::assign_kmer_order(uint16_t* kord, const char* seq, int k) {  // assumes a clean sequence (just 1s, 2s, 3s, 4s)

        size_t len = strlen(seq);
        if (len <= 0) throw std::runtime_error("Unexpected sequence length.");
        if (k >= len || k < 1 || k > 8) throw std::runtime_error("Invalid k-mer size.");

        size_t klen = len - k + 1; // number of k-mers in this sequence
        size_t n_kmers = (1 << (2 * k));  // 4^k k-mers
        if (kord == nullptr) throw std::runtime_error("Memory allocation failed.");
        for (int i = 0; i < klen; i++) {
            kord[i] = 0;
        }

        for (int i = 0; i < klen; i++) {

            size_t kmer = 0;
            for (int j = i; j < i + k; j++) {

                int nti = ((int)seq[j]) - 1; // change 1s, 2s, 3s, 4s, to 0/1/2/3
                if (nti != 0 && nti != 1 && nti != 2 && nti != 3) {
                    throw std::runtime_error("Unexpected nucleotide.");
                }

                kmer = 4 * kmer + nti;

            }

            // Make sure k-mer index is valid. This doesn't solve the N's/-'s
            // issue though, as the "length" of the string (# of k-mers) needs
            // to also reflect the reduction from the N's/-'s
            if (kmer == 999999) {
                ;
            } else if (kmer >= n_kmers) { // maybe check if above uin16_t max? (defined in header <stdint.h>, UINT16_MAX)
                throw std::runtime_error("Kmer index out of range.");
            } else { // valid k-mer
                kord[i] = kmer;
            }

        }

    }

    double Dada2Utility::kmer_dist(uint16_t* kv1, int len1, uint16_t* kv2, int len2, int k) {

        int n_kmer = 1 << (2 * k); // 4^k k-mers
        uint16_t dotsum = 0;

        for (int i = 0; i < n_kmer; i++) {
            dotsum += (kv1[i] < kv2[i] ? kv1[i] : kv2[i]);
        }

        double dot = ((double)dotsum) / ((len1 < len2 ? len1 : len2) - k + 1.0);

        return (1.0 - dot);

    }

    double Dada2Utility::kord_dist(uint16_t* kord1, int len1, uint16_t* kord2, int len2, int k) {

        uint16_t dotsum = 0;

        if (len1 != len2 || kord1 == nullptr || kord2 == nullptr) { // exit if different lengths
            return(-1.0);
        }
        size_t klen = len1 - k + 1;

        for (int i = 0; i < klen; i++) {
            dotsum += (kord1[i] == kord2[i]);
        }

        double dot = ((double)dotsum) / ((len1 < len2 ? len1 : len2) - k + 1.0);

        return (1.0 - dot);

    }


    char** Dada2Utility::nwalign_gapless(const char* s1, size_t len1, const char* s2, size_t len2) {

        size_t len_al = (len1 > len2) ? len1 : len2;

        // allocate memory to alignment strings
        char** al = new char*[2];
        al[0] = new char[len_al + 1];
        al[1] = new char[len_al + 1];

        // copy strings into the alignment strings
        for (int i = 0; i < len_al; i++) {

            al[0][i] = i < len1 ? s1[i] : '-';
            al[1][i] = i < len2 ? s2[i] : '-';

        }
        al[0][len_al] = '\0';
        al[1][len_al] = '\0';

        return al;

    }

    char** Dada2Utility::nwalign_endsfree(const char* s1, size_t len1, const char* s2, size_t len2, int score[4][4], int gap_p, int band) {

        unsigned int nrow = len1 + 1;
        unsigned int ncol = len2 + 1;
        int* d = new int[nrow * ncol];
        int* p = new int[nrow * ncol];

        // fill out left columns of d, p
        for (int i = 0; i <= len1; i++) {

            d[i * ncol] = 0; // ends-free gap
            p[i * ncol] = 3;

        }

        // fill out top rows of d, p
        for (int j = 0; j <= len2; j++) {

            d[j] = 0; // ends-free gap
            p[j] = 2;

        }

        // calculate left/right-bands in case of different lengths
        int lband, rband;
        if (len2 > len1) {

            lband = band;
            rband = band + len2 - len1;

        } else if (len1 > len2) {

            lband = band + len1 - len2;
            rband = band;

        } else {

            lband = band;
            rband = band;

        }

        // fill out band boundaries of d
        if (band >= 0 && (band < len1 || band < len2)) {

            for (int i = 0; i <= len1; i++) {

                if (i - lband - 1 >= 0) {
                    d[i * ncol + i - lband - 1] = -9999;
                }
                if (i + rband + 1 <= len2) {
                    d[i * ncol + i + rband + 1] = -9999;
                }

            }

        }

        // fill out the body of the DP matrix
        for (int i = 1; i <= len1; i++) {

            int l, r;

            if (band >= 0) {

                l = i - lband;
                if (l < 1) {
                    l = 1;
                }
                r = i + rband;
                if (r > len2) {
                    r = len2;
                }

            } else {

                l = 1;
                r = len2;

            }

            for (int j = l; j <= r; j++) {

                int left, up;

                // score for the left move
                if (i == len1) {
                    left = d[i * ncol + j - 1]; // ends-free gap
                } else {
                    left = d[i * ncol + j - 1] + gap_p;
                }

                // score for the up move
                if (j == len2) {
                    up = d[(i - 1) * ncol + j]; // ends-free gap
                } else {
                    up = d[(i - 1) * ncol + j] + gap_p;
                }

                // score for the diagonal move
                int diag = d[(i - 1) * ncol + j - 1] + score[s1[i - 1] - 1][s2[j - 1] - 1];

                // break ties and fill in d, p
                if (up >= diag && up >= left) {

                    d[i * ncol + j] = up;
                    p[i * ncol + j] = 3;

                } else if (left >= diag) {

                    d[i * ncol + j] = left;
                    p[i * ncol + j] = 2;

                } else {

                    d[i * ncol + j] = diag;
                    p[i * ncol + j] = 1;

                }

            }

        }

        char* al0 = new char[len1 + len2 + 1];
        char* al1 = new char[len1 + len2 + 1];

        // trace back over p to form the alignment
        size_t len_al = 0;
        int i = len1;
        int j = len2;

        while (i > 0 || j > 0) {

            switch (p[i * ncol + j]) {
                case 1:
                    al0[len_al] = s1[--i];
                    al1[len_al] = s2[--j];
                    break;
                case 2:
                    al0[len_al] = '-';
                    al1[len_al] = s2[--j];
                    break;
                case 3:
                    al0[len_al] = s1[--i];
                    al1[len_al] = '-';
                    break;
                default:
                    throw std::runtime_error("N-W Align out of range.");
            }
            len_al++;

        }
        al0[len_al] = '\0';
        al1[len_al] = '\0';


        // allocate memory to alignment strings
        char** al = new char*[2];
        al[0] = new char[len_al + 1];
        al[1] = new char[len_al + 1];

        // reverse the alignment strings (since traced backwards)
        for (i = 0; i < len_al; i++) {

            al[0][i] = al0[len_al - i - 1];
            al[1][i] = al1[len_al - i - 1];

        }
        al[0][len_al] = '\0';
        al[1][len_al] = '\0';

        // free allocated memory
        delete[] d;
        delete[] p;
        delete[] al0;
        delete[] al1;

        return al;

    }

    char** Dada2Utility::nwalign_endsfree_homo(const char* s1, size_t len1, const char* s2, size_t len2, int score[4][4], int gap_p, int homo_gap_p, int band) {

        // find locations where s1 has homopolymer and put 1s in homo1
        unsigned char* homo1 = new unsigned char[len1];
        unsigned char* homo2 = new unsigned char[len2];

        for (int i = 0, j = 0; j < len1; j++) {

            if (j == len1 - 1 || s1[j] != s1[j + 1]) {

                for (int k = i; k <= j; k++) {

                    if (j - i >= 2) {// min homopolymer length = 3
                        homo1[k] = 1;
                    } else {
                        homo1[k] = 0;
                    }

                }
                i = j + 1;

            }

        }

        // find locations where s2 has homopolymer and put 1s in homo2
        for (int i = 0, j = 0; j < len2; j++) {

            if (j == len2 - 1 || s2[j] != s2[j + 1]) {

                for (int k = i; k <= j; k++) {

                    if (j - i >= 2) { // min homopolymer length = 3
                        homo2[k] = 1;
                    } else {
                        homo2[k] = 0;
                    }

                }
                i = j + 1;

            }

        }

        unsigned int nrow = len1 + 1;
        unsigned int ncol = len2 + 1;
        int* d = new int[nrow * ncol];
        int* p = new int[nrow * ncol];

        // fill out left columns of d, p
        for (int i = 0; i <= len1; i++) {

            d[i * ncol] = 0; // ends-free gap
            p[i * ncol] = 3;

        }

        // fill out top rows of d, p
        for (int j = 0; j <= len2; j++) {
            d[j] = 0; // ends-free gap
            p[j] = 2;
        }

        // calculate left/right-bands in case of different lengths
        int lband, rband;
        if (len2 > len1) {

            lband = band;
            rband = band + len2 - len1;

        } else if (len1 > len2) {

            lband = band + len1 - len2;
            rband = band;

        } else {

            lband = band;
            rband = band;

        }

        // fill out band boundaries of d
        if (band >= 0 && (band < len1 || band < len2)) {

            for (int i = 0; i <= len1; i++) {

                if (i - lband - 1 >= 0) {
                    d[i * ncol + i - lband - 1] = -9999;
                }
                if (i + rband + 1 <= len2) {
                    d[i * ncol + i + rband + 1] = -9999;
                }

            }

        }

        // fill out the body of the DP matrix
        for (int i = 1; i <= len1; i++) {

            int l, r;

            if (band >= 0) {

                l = i - lband;
                if (l < 1) {
                    l = 1;
                }
                r = i + rband;
                if (r > len2) {
                    r = len2;
                }

            } else {

                l = 1;
                r = len2;

            }

            for (int j = l; j <= r; j++) {

                int left, up;

                // score for the left move
                if (i == len1) {
                    left = d[i * ncol + j - 1]; // ends-free gap
                } else if (homo2[j - 1]) {
                    left = d[i * ncol + j - 1] + homo_gap_p; // homopolymer gap
                } else {
                    left = d[i * ncol + j - 1] + gap_p;
                }

                // score for the up move
                if (j == len2) {
                    up = d[(i - 1) * ncol + j]; // ends-free gap
                } else if (homo1[i-1]) {
                    up = d[(i - 1) * ncol + j] + homo_gap_p; // homopolymer gap
                } else {
                    up = d[(i - 1) * ncol + j] + gap_p;
                }

                // score for the diagonal move
                int diag = d[(i - 1) * ncol + j - 1] + score[s1[i - 1] - 1][s2[j - 1] - 1];

                // break ties and fill in d, p
                if (up >= diag && up >= left) {

                    d[i * ncol + j] = up;
                    p[i * ncol + j] = 3;
                } else if (left >= diag) {

                    d[i * ncol + j] = left;
                    p[i * ncol + j] = 2;
                } else {

                    d[i * ncol + j] = diag;
                    p[i * ncol + j] = 1;

                }
            }
        }

        char* al0 = new char[len1 + len2 + 1];
        char* al1 = new char[len1 + len2 + 1];

        // trace back over p to form the alignment
        size_t len_al = 0;
        int i = len1;
        int j = len2;

        while (i > 0 || j > 0) {

            switch (p[i * ncol + j]) {
                case 1:
                    al0[len_al] = s1[--i];
                    al1[len_al] = s2[--j];
                    break;
                case 2:
                    al0[len_al] = '-';
                    al1[len_al] = s2[--j];
                    break;
                case 3:
                    al0[len_al] = s1[--i];
                    al1[len_al] = '-';
                    break;
                default:
                    throw std::runtime_error("N-W Align out of range.");
            }
            len_al++;

        }
        al0[len_al] = '\0';
        al1[len_al] = '\0';


        // allocate memory to alignment strings
        char**al = new char*[2];
        al[0] = new char[len_al + 1];
        al[1] = new char[len_al + 1];

        // reverse the alignment strings (since traced backwards)
        for (int ii = 0; ii < len_al; ii++) {

            al[0][ii] = al0[len_al - ii - 1];
            al[1][ii] = al1[len_al - ii - 1];

        }
        al[0][len_al] = '\0';
        al[1][len_al] = '\0';

        // free allocated memory
        delete[] d;
        delete[] p;
        delete[] homo1;
        delete[] homo2;
        delete[] al0;
        delete[] al1;

        return al;

    }


    double Dada2Utility::calc_p_a(int reads, double exp_reads, bool prior) {

        // calculate p-value from Poisson cdf
        numSeqs_t n_repeats = reads - 1; // -1 since strict > being calculated, and want to include the observed count
        double pval = 1 - stats::ppois(n_repeats, exp_reads, false); // lower.tail = false: P(X > x) --> here: 1 - P(X <= x)

        if (!prior) {

            // calculate norm (since conditioning on sequence being present)
            double norm = (1.0 - exp(-exp_reads));
            if (norm < tail_approx_cutoff) {
                norm = exp_reads - 0.5 * exp_reads * exp_reads;
                // assumption: TAIL_APPROX_CUTOFF is small enough to terminate taylor expansion at 2nd order
            }
            pval = pval / norm;

        }

        return pval;

    }


    Dada2Utility::RawSequence::RawSequence() {

        seq = nullptr;
        qual = nullptr;
        prior = false;
        kmer = nullptr;
        kmer8 = nullptr;
        kord = nullptr;
        length = 0;
        abundance = 0;
        index = 0;
        p = 0.0;
        max_exp = -999.0;
        lock = false;

    }

    Dada2Utility::RawSequence::RawSequence(const char* ampl_seq, lenSeqs_t ampl_len, uint8_t* ampl_quals, unsigned int ampl_abundance, bool is_prior) {

        // assign sequence and associated properties
        length = ampl_len;
        abundance = ampl_abundance;
        seq = new char[length + 1];
        strcpy(seq, ampl_seq);
        qual = new uint8_t[length];
        memcpy(qual, ampl_quals, length * sizeof(uint8_t));

        prior = is_prior;
        p = 0.0;
        max_exp = -999.0;
        lock = false;

        kmer = nullptr;
        kmer8 = nullptr;
        kord = nullptr;
        index = 0;

    }

    Dada2Utility::RawSequence::~RawSequence() {

        delete[] qual;
        delete[] seq;

    }

    Dada2Utility::RawSequence::RawSequence(const RawSequence& other) { // copy constructor

        length = other.length;
        abundance = other.abundance;
        prior = other.prior;
        index = other.index;
        p = other.p;
        max_exp = other.max_exp;
        comp = other.comp;
        lock = other.lock;

        seq = new char[length + 1];
        strcpy(seq, other.seq);
        qual = new uint8_t[length];
        memcpy(qual, other.qual, length * sizeof(uint8_t));

        size_t n_kmer = 1 << (2 * Dada2Utility::k_mer_size);
        if (other.kmer == nullptr) {
            kmer = nullptr;
        } else {
            kmer = new uint16_t[n_kmer];
            memcpy(kmer, other.kmer, n_kmer * sizeof(uint16_t));
        }
        if (other.kmer8 == nullptr) {
            kmer8 = nullptr;
        } else {
            kmer8 = new uint8_t[n_kmer];
            memcpy(kmer8, other.kmer8, n_kmer * sizeof(uint8_t));
        }
        if (other.kord == nullptr) {
            kord = nullptr;
        } else {
            kord = new uint16_t[length];
            memcpy(kord, other.kord, length * sizeof(uint16_t));
        }

    }

    Dada2Utility::RawSequence::RawSequence(RawSequence&& other) noexcept { // move constructor

        length = other.length;
        abundance = other.abundance;
        prior = other.prior;
        index = other.index;
        p = other.p;
        max_exp = other.max_exp;
        comp = other.comp;
        lock = other.lock;

        seq = other.seq; other.seq = nullptr;
        qual = other.qual; other.qual = nullptr;

        kmer = other.kmer; other.kmer = nullptr;
        kmer8 = other.kmer8; other.kmer8 = nullptr;
        kord = other.kord; other.kord = nullptr;

    }

    Dada2Utility::RawSequence& Dada2Utility::RawSequence::operator=(const RawSequence& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] seq;
        delete[] qual;
        delete[] kmer;
        delete[] kmer8;
        delete[] kord;

        // copy new resources
        length = other.length;
        abundance = other.abundance;
        prior = other.prior;
        index = other.index;
        p = other.p;
        max_exp = other.max_exp;
        comp = other.comp;
        lock = other.lock;

        seq = new char[length + 1];
        strcpy(seq, other.seq);
        qual = new uint8_t[length];
        memcpy(qual, other.qual, length * sizeof(uint8_t));

        size_t n_kmer = 1 << (2 * Dada2Utility::k_mer_size);
        if (other.kmer == nullptr) {
            kmer = nullptr;
        } else {
            kmer = new uint16_t[n_kmer];
            memcpy(kmer, other.kmer, n_kmer * sizeof(uint16_t));
        }
        if (other.kmer8 == nullptr) {
            kmer8 = nullptr;
        } else {
            kmer8 = new uint8_t[n_kmer];
            memcpy(kmer8, other.kmer8, n_kmer * sizeof(uint8_t));
        }
        if (other.kord == nullptr) {
            kord = nullptr;
        } else {
            kord = new uint16_t[length];
            memcpy(kord, other.kord, length * sizeof(uint16_t));
        }

        return *this;

    }

    Dada2Utility::RawSequence& Dada2Utility::RawSequence::RawSequence::operator=(RawSequence&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] seq;
        delete[] qual;
        delete[] kmer;
        delete[] kmer8;
        delete[] kord;

        // copy / transfer new resources
        length = other.length;
        abundance = other.abundance;
        prior = other.prior;
        index = other.index;
        p = other.p;
        max_exp = other.max_exp;
        comp = other.comp;
        lock = other.lock;

        seq = other.seq; other.seq = nullptr;
        qual = other.qual; other.qual = nullptr;

        kmer = other.kmer; other.kmer = nullptr;
        kmer8 = other.kmer8; other.kmer8 = nullptr;
        kord = other.kord; other.kord = nullptr;

        return *this;

    }


    char** Dada2Utility::Sub::raw_align(RawSequence* raw1, RawSequence* raw2, int match, int mismatch, int gap_p, int homo_gap_p,
                                        bool use_k_mers, double k_dist_cutoff, int band, bool gapless) {

        // k-mer screen
        double k_dist = 0.0;
        if (use_k_mers) {
            k_dist = kmer_dist(raw1->kmer, raw1->length, raw2->kmer, raw2->length, k_mer_size); // implicit vectorisation
        }

        // gapless screen (using k-mers)
        double ko_dist = -1.0; // needs to be different than kdist for fall-back when use_k_mers=false
        if (use_k_mers && gapless) {
            ko_dist = kord_dist(raw1->kord, raw1->length, raw2->kord, raw2->length, k_mer_size);
        }

        // make deprecated score matrix
        int score[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                score[i][j] = i==j ? match : mismatch;
            }
        }

        char** al;
        if (use_k_mers && k_dist > k_dist_cutoff) {
            al = nullptr;
        } else if (band == 0 || (gapless && ko_dist == k_dist)) {
            al = nwalign_gapless(raw1->seq, raw1->length, raw2->seq, raw2->length);
        } else if (homo_gap_p != gap_p && homo_gap_p <= 0) {
            al = nwalign_endsfree_homo(raw1->seq, raw1->length, raw2->seq, raw2->length, score, gap_p, homo_gap_p, band); // uses old score-matrix format
        } else {
            al = nwalign_endsfree(raw1->seq, raw1->length, raw2->seq, raw2->length, score, gap_p, band); // uses old score-matrix format
        }

        return al;

    }

    Dada2Utility::Sub::Sub(RawSequence* raw0, RawSequence* raw1, int match, int mismatch, int gap_p, int homo_gap_p, bool use_kmers,
                           double kdist_cutoff, int band, bool gapless) {

        char** al = raw_align(raw0, raw1, match, mismatch, gap_p, homo_gap_p, use_kmers, kdist_cutoff, band, gapless);

        len0 = 0;
        nsubs = 0;
        map = nullptr;
        pos = nullptr;
        nt0 = nullptr;
        nt1 = nullptr;

        // null alignment (outside k-mer threshold) -> null sub; otherwise: process the alignment
        if (al) {

            // traverse alignment and find length of seq0 and nsubs for memory allocation
            char* al0 = al[0]; // dummy pointers to the sequences in the alignment
            char* al1 = al[1]; // dummy pointers to the sequences in the alignment
            int align_length = strlen(al[0]);
            for (int i = 0; i < align_length; i++) {

                bool is_nt0 = ((al0[i] == 1) || (al0[i] == 2) || (al0[i] == 3) || (al0[i] == 4) || (al0[i] == 5)); // A/C/G/T/N (non-gap) in seq0
                bool is_nt1 = ((al1[i] == 1) || (al1[i] == 2) || (al1[i] == 3) || (al1[i] == 4) || (al1[i] == 5)); // A/C/G/T/N (non-gap) in seq1

                if (is_nt0) {
                    len0++;
                }
                if (is_nt0 && is_nt1) { // possible sub
                    if ((al0[i] != al1[i]) && (al0[i] != 5) && (al1[i] != 5)) { // Ns don't make subs
                        nsubs++;
                    }
                }

            }

            map = new uint16_t[len0];
            pos = new uint16_t[nsubs];
            nt0 = new char[nsubs];
            nt1 = new char[nsubs];
            nsubs = 0;

            // traverse the alignment and record substitutions
            int i0 = -1;
            int i1 = -1;
            al0 = al[0];
            al1 = al[1];
            for (int i = 0; i < align_length; i++) {

                bool is_nt0 = ((al0[i] == 1) || (al0[i] == 2) || (al0[i] == 3) || (al0[i] == 4) || (al0[i] == 5)); // A/C/G/T/N (non-gap) in seq0
                bool is_nt1 = ((al1[i] == 1) || (al1[i] == 2) || (al1[i] == 3) || (al1[i] == 4) || (al1[i] == 5)); // A/C/G/T/N (non-gap) in seq1
                if (is_nt0) {
                    i0++;
                }
                if (is_nt1) {
                    i1++;
                }

                // record to map
                if (is_nt0) {

                    if (is_nt1) {
                        map[i0] = i1; // assigning a signed into to a uint16_t...
                    } else {
                        map[i0] = gap_glyph; // indicates gap
                    }

                }

                if (is_nt0 && is_nt1) { // possible sub

                    if ((al0[i] != al1[i]) && (al0[i] != 5) && (al1[i] != 5)) { // Ns don't make subs

                        pos[nsubs] = i0;
                        nt0[nsubs] = al0[i];
                        nt1[nsubs] = al1[i];
                        nsubs++;

                    }

                }

            }

            delete[] al[0];
            delete[] al[1];
            delete[] al;

        }

    }

    Dada2Utility::Sub::~Sub() {

        delete[] nt1;
        delete[] nt0;
        delete[] pos;
        delete[] map;

    }

    double Dada2Utility::compute_lambda(const char* seq, size_t seq_len, const uint8_t* qual, Sub* sub, std::vector<std::vector<double>>& error_matrix, int* tvec, unsigned int* qind) {

        if ((sub == nullptr) || (sub->nsubs == 0 && sub->len0 == 0)) { // null Sub / outside k-mer threshold
            return 0.0;
        }

        // make vector that indexes as integers the transitions at each position in seq1
        // index is 0: A->A, 1: A->C, ..., 4: C->A, ...
        for (int pos1 = 0; pos1 < seq_len; pos1++) {

            int nti1 = ((int)seq[pos1]) - 1;
            if (nti1 == 0 || nti1 == 1 || nti1 == 2 || nti1 == 3) {
                tvec[pos1] = nti1 * 4 + nti1;
            } else {
                throw std::runtime_error("Non-ACGT sequences in compute_lambda.");
            }

            // turn quality into the index in the array
            qind[pos1] = qual[pos1];

        }

        // now fix the positions where substitutions occurred
        for (int s = 0; s < sub->nsubs; s++) {

            int pos0 = sub->pos[s];
            if (pos0 < 0 || pos0 >= sub->len0) throw std::runtime_error("CL: Bad pos0: " + std::to_string(pos0) + "(len0=" + std::to_string(sub->len0) + ").");
            int pos1 = sub->map[sub->pos[s]];
            if (pos1 < 0 || pos1 >= seq_len) throw std::runtime_error("CL: Bad pos1: " + std::to_string(pos1) + " (len1=" + std::to_string(seq_len) + ").");

            int nti0 = ((int)sub->nt0[s]) - 1;
            int nti1 = ((int)sub->nt1[s]) - 1;
            tvec[pos1] = nti0 * 4 + nti1;

        }

        // calculate lambda
        double lambda = 1.0;
        for (int pos1 = 0; pos1 < seq_len; pos1++) {
            lambda = lambda * error_matrix.at(tvec[pos1]).at(qind[pos1]);
        }

        if (lambda < 0 || lambda > 1) throw std::runtime_error("Bad lambda.");

        return lambda;

    }

    bool Dada2Utility::fill_error_matrix(std::vector<std::vector<double>>& error_matrix, const std::string& file, short num_scores) {

        if (file.empty()) return false;

        std::ifstream in_stream(file);
        if (!in_stream.good()) {

            std::cerr << " -- WARNING: File '" << file << "' not opened correctly. Error matrix could not be read from it." << std::endl;
            return false;

        }

        std::string line;
        bool was_extended = false;
        while (std::getline(in_stream, line).good()) {

            error_matrix.emplace_back();
            auto& row = error_matrix.back();

            // fill (beginning of) matrix row with values from the file
            while (line.length() > 0) {

                std::string val = line.substr(0, line.find(' '));
                row.push_back(std::stod(val));
                line.erase(0, val.length() + 1);

            }

            was_extended |= (row.size() < num_scores);

            // extend matrix row to the number of needed scores by repeating the last value
            while (row.size() < num_scores) {
                row.push_back(row.back());
            }

        }

        if (was_extended) {
            std::cerr << " -- WARNING: The provided error matrix did not have enough columns and was extended by repeating the last provided value per line." << std::endl;
        }

        return true;

    }



}
