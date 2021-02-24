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

#ifndef GEFAST_UTILITY_HPP
#define GEFAST_UTILITY_HPP

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "Options.hpp"

namespace GeFaST {

    // version number of current implementation of GeFaST
    const std::string VERSION = "2.0.1";

    // output separators
    const char SEP_INTERNALS = '\t'; // between columns in output_internal
    const char SEP_OTUS = ' '; // between amplicon identifiers in output_otus
    const char SEP_STATISTICS = '\t'; // between columns in output_statistics
    const char SEP_UCLUST = '\t'; // between columns in output_uclust
    const char SEP_MOTHUR = ','; // for mothur-compatible output
    const char SEP_MOTHUR_OTU = '\t'; // for mothur-compatible output

    // type for everything related to counts of amplicons / sequences
    typedef unsigned long numSeqs_t;

    // type for everything related to the length of the sequence of an amplicon
    typedef unsigned long lenSeqs_t;

    // type for everything related to the distance between two amplicons
    typedef float dist_t;



    /*
     * Print general information on the tool as a header.
     */
    void print_information();


    /*
     * Print a message briefly summarizing the command-line options.
     */
    void print_help();


    /*
     * Read file expected contain one file path per line.
     * Empty lines and comment lines (starting with #) are allowed.
     */
    std::vector<std::string> read_list_file(const std::string& list_file);


    /*
     * Check whether the given sequence satisfies the given filter criteria.
     */
    bool check_sequence(const std::string& seq, const std::string& alphabet, numSeqs_t seq_abund, lenSeqs_t min_len,
            lenSeqs_t max_len, numSeqs_t min_abund, numSeqs_t max_abund);


    /*
     * Set sequence in upper case.
     */
    void upper_case(std::string& s);


    /*
     * Compute gcd of two non-negative integers.
     */
    unsigned long long gcd(unsigned long long a, unsigned long long b);


    /*
     * Representation of a header / defline of a sequence entry.
     */
    struct Defline {

        std::string id = ""; // string identifier of amplicon
        std::string information = ""; // additional information (2nd part of defline)
        numSeqs_t abundance = 0; // abundance of amplicon

        /*
         * Construct 'empty' instance with no empty id, empty info and abundance of 0.
         */
        Defline() = default;

        /*
         * Construct 'filled' instance with the given identifier i, the additional info and abundance a.
         */
        Defline(const std::string& i, const std::string& info, const numSeqs_t a);

        /*
         * Split description line into actual header, abundance value and additional information (if any).
         * Suitable for FASTA and FASTQ files.
         */
        static Defline parse_description_line(const std::string& defline, const std::string& sep);

    };


    /*
     * Class for recording statistics on the input.
     *
     * Records the following information:
     *  - total number of amplicons
     *  - number of amplicons per sequence length
     *  - total length of headers (identifiers) of all amplicons
     *  - total length of headers (identifiers) of amplicons per sequence length
     *  - total length of sequences of all amplicons
     *
     *  Offers methods for recording and unrecord amplicons provided in different forms.
     *  Only considers the lengths of headers and sequences, not their content.
     */
    template<typename L = lenSeqs_t, typename N = numSeqs_t>
    class DataStatistics {

    public:
        /* === Add / remove information === */

        void record(const std::string header, const std::string seq) {

            total_num++;
            total_length_headers += header.length();
            total_length_seqs += seq.length();
            num_per_length[seq.length()]++;
            total_length_headers_per_length[seq.length()] += header.length();

        }

        void record(const L len_header, const L len_seq) {

            total_num++;
            total_length_headers += len_header;
            total_length_seqs += len_seq;
            num_per_length[len_seq]++;
            total_length_headers_per_length[len_seq] += len_header;

        }

        void record(const char* header, const L len_header, const char* seq, const L len_seq) {

            total_num++;
            total_length_headers += len_header;
            total_length_seqs += len_seq;
            num_per_length[len_seq]++;
            total_length_headers_per_length[len_seq] += len_header;

        }

        void unrecord(const std::string header, const std::string seq) {

            total_num--;
            total_length_headers -= header.length();
            total_length_seqs -= seq.length();
            num_per_length[seq.length()]--;
            total_length_headers_per_length[seq.length()] -= header.length();

        }

        void unrecord(const L len_header, const L len_seq) {

            total_num--;
            total_length_headers -= len_header;
            total_length_seqs -= len_seq;
            num_per_length[len_seq]--;
            total_length_headers_per_length[len_seq] -= len_header;

        }

        void unrecord(const char* header, const L len_header, const char* seq, const L len_seq) {

            total_num--;
            total_length_headers -= len_header;
            total_length_seqs -= len_seq;
            num_per_length[len_seq]--;
            total_length_headers_per_length[len_seq] -= len_header;

        }


        /* === Retrieve information === */

        unsigned long long get_total_length_headers() const {
            return total_length_headers;
        }

        unsigned long long get_total_length_sequences() const {
            return total_length_seqs;
        }

        N get_total_num() const {
            return total_num;
        }

        N get_num_per_length(const L len) const {

            auto iter = num_per_length.find(len);
            return (iter != num_per_length.end()) ? iter->second : 0;

        }

        N get_total_length_headers_per_length(const L len) const {

            auto iter = total_length_headers_per_length.find(len);
            return (iter != total_length_headers_per_length.end()) ? iter->second : 0;

        }

        L get_min_length() const {
            return (num_per_length.empty()) ? 0 : num_per_length.cbegin()->first;
        }

        L get_max_length() const {
            return (num_per_length.empty()) ? 0 : num_per_length.crbegin()->first;
        }

        // Return the lengths in ascending order.
        std::vector<L> get_all_lengths() const {

            std::vector<L> lengths;
            for (auto iter = num_per_length.cbegin(); iter != num_per_length.cend(); iter++) {
                lengths.push_back(iter->first);
            }

            return lengths;

        }

    protected:
        N total_num = 0; // total number of amplicons
        unsigned long long total_length_headers = 0; // total length of headers (identifiers) of all amplicons
        unsigned long long total_length_seqs = 0; // total length of sequences of all amplicons
        std::map<L, N> num_per_length; // number of amplicons per observed sequence length
        std::map<L, unsigned long long> total_length_headers_per_length; // total length of headers (identifiers) of amplicons per observed sequence length

    };


    /*
     * Class for converting quality scores to error probabilities and vice versa.
     *
     * Currently, 4 encodings are supported:
     *  - Sanger (phred+33, 0 to 40)
     *  - Illumina 1.3+ (phred+64, 0 to 40)
     *  - Illumina 1.5+ (phred+64, 3 to 41)
     *  - Illumina 1.8+ (phred+33, 0 to 41)
     */
    template<typename T = float>
    class QualityEncoding {

    public:
        explicit QualityEncoding(const QualityEncodingOption opt = QE_SANGER) {

            int from, to;

            switch (opt) {

                case QE_SANGER:
                    offset_ = 33;
                    from = 0;
                    to = 40;
                    break;

                case QE_ILLUMINA13:
                    offset_ = 64;
                    from = 0;
                    to = 40;
                    break;

                case QE_ILLUMINA15:
                    offset_ = 64;
                    from = 3;
                    to = 41;
                    break;

                case QE_ILLUMINA18:
                    offset_ = 33;
                    from = 0;
                    to = 41;
                    break;

                case QE_UNKNOWN: // fallthrough to default for error handling
                default:
                    std::cerr << " -- WARNING: Unknown quality-encoding name. Resorting to sanger." << std::endl;
                    offset_ = 33;
                    from = 0;
                    to = 40;

            }

            for (char i = from; i <= to; i++) {
                score_to_error_prob_[i + offset_] = quality_to_probability(i);
            }

        }

        /*
         * Convert the quality score to the error probability.
         */
        inline T quality_to_probability(int quality) const {
            return powf(10, (quality / (T)-10.0));
        }

        /*
         * Convert the error probability to the quality score.
         */
        inline int probability_to_quality(T probability) const {
            return -10 * log10(probability);
        }

        /*
         * Convert the encoded quality score (contains offset for visibility as ASCII character)
         * to the error probability.
         */
        inline T encoded_quality_to_probability(char quality) const {
            //return pow(10, ((quality - offset_) / (T) -10.0));
            return score_to_error_prob_.at(quality);
        }

        /*
         * Convert the error probability to the quality score and add an offset to
         * obtain a visible ASCII character (encoded quality score).
         */
        inline char probability_to_encoded_quality(T probability) const {
            return std::round(-10 * log10(probability) + offset_);
        }

        /*
         * List all quality-score values that are accepted in the quality encoding.
         */
        std::vector<char> get_accepted_scores() const {

            std::vector<char> scores;
            for (auto iter = score_to_error_prob_.begin(); iter != score_to_error_prob_.end(); iter++) {
                scores.push_back(iter->first);
            }

            return scores;

        }

        /*
         * Return the offset value of the quality encoding.
         */
        char get_offset() const {
            return offset_;
        }

    private:
        std::map<char, T> score_to_error_prob_; // mapping between quality score and error probability
        char offset_; // encoding-specific offset

    };


    /*
     * Class encapsulating further information that might be necessary for some modes, e.g. quality scores.
     * These information are handed over to the AmpliconStorage together with the other amplicon information
     * during the preprocessing.
     */
    struct ExtraInfo {
        // empty base class in case there are no information
    };

    /*
     * Further information on amplicon for modes requiring quality scores.
     *
     * Provides the quality scores as a string of encoded ASCII characters.
     */
    struct ExtraInfoQuality : public ExtraInfo {

        ExtraInfoQuality() = default;

        std::string quality_scores = ""; // encoded quality scores
    };


    /*
     * Pair of pointers (first, second) describing the (sub)string corresponding to [first, last).
     */
    typedef std::pair<const char*, const char*> StringIteratorPair;

    // hash function for StringIteratorPair
    struct hashStringIteratorPair {

        hashStringIteratorPair() : hash(std::hash<std::string>()) {
            // nothing else to do
        }

        size_t operator()(const StringIteratorPair& p) const;

        std::hash<std::string> hash;
    };

    // comparison function for string equality
    struct equalStringIteratorPair {
        bool operator()(const StringIteratorPair& lhs, const StringIteratorPair& rhs) const;
    };

    // comparison function for lexicographical order
    struct lessStringIteratorPair {
        bool operator()(const StringIteratorPair& a, const StringIteratorPair& b) const;
    };


    // hash function for std::pair of integers
    // based on: https://www.boost.org/doc/libs/1_35_0/doc/html/boost/hash_combine_id241013.html
    struct hashStdPair {

        hashStdPair() : hash(std::hash<lenSeqs_t>()) {
            // nothing else to do
        }

        size_t operator()(const std::pair<lenSeqs_t, lenSeqs_t>& p) const;

        std::hash<lenSeqs_t> hash;

    };


    /*
     * Utility methods and variables used by the DADA2-inspired, consistency-based
     * clustering and refinement methods.
     *
     * Most of the methods are slight adaptations of eponymous methods in the DADA2 implementation:
     * "DADA2: High-resolution sample inference from Illumina amplicon data"
     * (Callahan et al., 2016, https://doi.org/10.1038/nmeth.3869)
     *
     * Implementation obtained from: https://github.com/benjjneb/dada2
     */
    struct Dada2Utility {


        static const int k_mer_size = 5; // length of k-mers considered in distance screen prior to alignment computation (if activated)
        static const int gap_glyph = 9999; // internal representation of gap symbol
        static constexpr double tail_approx_cutoff = 1e-7; // cutoff for tail approximation P(X >= 1) of Poisson distribution
        static const int raw_buffer = 50; // capacity increment (number of RawSequence instances) when allocated space is insufficient
        static const int clust_buffer = 50; // capacity increment (number of Cluster instances) when allocated space is insufficient
        static const int max_shuffle = 10; // maximum number of shuffle rounds after a new cluster has budded


        // Convert an input string (iseq) to numeric indices, which are stored in the char output array (oseq).
        // Needs memory at oseq equal to strlen of iseq, which is not checked currently.
        // Source: nt2int from misc.cpp
        static void nt2int(char* oseq, const char* iseq);

        // Convert sequence in index form back to nucleotides.
        // Source: int2nt from misc.cpp
        static void int2nt(char* oseq, const char* iseq);


        // Compute k-mer profile stored in 8-bit integers.
        // The counts of k-mers occurring more than 255 times are capped at 255.
        // Source: assign_kmer8 from kmers.cpp
        static void assign_kmer8(uint8_t* kvec8, const char* seq, int k);

        // Compute k-mer profile stored in 16-bit integers.
        // Source: assign_kmer from kmers.cpp
        static void assign_kmer(uint16_t* kvec, const char* seq, int k);

        // Compute the vector of the k-mers in order of appearance.
        // Source: assign_kmer_order from kmers.cpp
        static void assign_kmer_order(uint16_t* kord, const char* seq, int k);

        // Compute k-mer distance based on k-mer counts.
        // Source: kmer_dist from kmers.cpp
        static double kmer_dist(uint16_t* kv1, int len1, uint16_t* kv2, int len2, int k);

        // Compute the k-mer distance based on ordered k-mers. Return -1 (invalid) if the lengths differ.
        // Source: kord_dist from kmers.cpp
        static double kord_dist(uint16_t* kord1, int len1, uint16_t* kord2, int len2, int k);


        // Compute gapless alignment (gap symbols can only occur at the end of the shorter sequence).
        // Source: nwalign_gapless from nwalign_endsfree.cpp
        static char** nwalign_gapless(const char* s1, size_t len1, const char* s2, size_t len2);

        // Compute alignment with free gap symbols at the ends of the sequences.
        // Note: input sequence must end with string termination character '\0'.
        // Source: nwalign_endsfree from nwalign_endsfree.cpp
        static char** nwalign_endsfree(const char* s1, size_t len1, const char* s2, size_t len2, int score[4][4], int gap_p, int band);

        // Compute alignment with free gap symbols at the ends of the sequence
        // and a special alignment gap penalty within homopolymer regions.
        // Note: input sequence must end with string termination character '\0'.
        // Source: nwalign_endsfree_homo from nwalign_endsfree.cpp
        static char** nwalign_endsfree_homo(const char* s1, size_t len1, const char* s2, size_t len2, int score[4][4], int gap_p, int homo_gap_p, int band);


        // Calculate abundance p-value for given abundance and expected abundance.
        // The p-value is conditional on the sequence being present, unless prior evidence is true.
        // Source: calc_pA from pval.cpp
        static double calc_p_a(int reads, double exp_reads, bool prior);


        // Summary of comparison between a cluster and an amplicon.
        // Source: Comparison from dada.h
        struct Comparison {

            unsigned int i = 0;         // index of cluster
            unsigned int index = 0;     // index of raw sequence
            double lambda = 0.0;        // error model lambda
            unsigned int hamming = 0;   // number of substitutions

        };

        // Representation of an amplicon with additional information.
        // Source: Raw from dada.h
        struct RawSequence {

            char* seq;              // the sequence, stored as C-string with A = 1, C = 2, G = 3, T = 4
            uint8_t* qual;          // the rounded average qualities at each position for this unique
            bool prior;             // there are (not) prior reasons to expect this sequence to exist
            uint16_t* kmer;         // the 16-bit k-mer vector (counts) of this sequence
            uint8_t* kmer8;         // the 8-bit k-mer vector (counts) of this sequence
            uint16_t* kord;         // the k-mers of this sequence in their order of appearance
            unsigned int length;    // the length of the sequence
            unsigned int abundance; // abundance (number of reads) of this (unique) sequence
            unsigned int index;     // the index of this raw in a collection of raws (e.g. ClusterCollection::raws)
            double p;               // abundance p-value relative to the current cluster
            double max_exp;         // minimum expected number of reads, maximised over examined clusters
            Comparison comp;        // the comparison between this Raw and...
            bool lock;              // "locks" the raw to its current cluster in greedy-mode


            // empty constructor for null initialisation
            RawSequence();

            // Basic constructor. Does not fill the k-mer vectors.
            // Source: raw_new from containers.cpp
            RawSequence(const char* ampl_seq, lenSeqs_t ampl_len, uint8_t* ampl_quals, unsigned int ampl_abundance, bool is_prior);

            ~RawSequence();

            RawSequence(const RawSequence& other); // copy constructor

            RawSequence(RawSequence&& other) noexcept; // move constructor

            RawSequence& operator=(const RawSequence& other); // copy assignment operator

            RawSequence& operator=(RawSequence&& other) noexcept; // move assignment operator

        };

        // Representation of a set of substitutions (position and identity) of one sequence in an alignment to another sequence.
        // Source: Sub from dada.h
        struct Sub {

            unsigned int nsubs; // number of substitutions
            unsigned int len0;  // length of the alignment's reference sequence
            uint16_t* map;      // map of the sequence position in the reference sequence to that in the aligned sequence
            uint16_t* pos;      // sequence position of the substitution: index in the reference sequence
            char* nt0;          // nucleotide in reference sequence
            char* nt1;          // different nucleotide in aligned sequence

            // Alignment computation (banded Needleman-Wunsch algorithm)
            // Source: raw_align from nwalign_endsfree.cpp
            char** raw_align(RawSequence* raw1, RawSequence* raw2, int match, int mismatch, int gap_p, int homo_gap_p,
                    bool use_k_mers, double k_dist_cutoff, int band, bool gapless);

            // Constructor computing the alignment and extracting the information on the substitutions.
            // Source: sub_new, al2subs from nwalign_endsfree.cpp
            Sub(RawSequence* raw0, RawSequence* raw1, int match, int mismatch, int gap_p, int homo_gap_p, bool use_kmers,
                double kdist_cutoff, int band, bool gapless);

            ~Sub();

        };

        // Calculate error rate (lambda) from a lookup table indexed by transition (row) and rounded quality (col)
        // using a sequence and a set of substitutions.
        // Source: compute_lambda from pval.cpp
        static double compute_lambda(const char* seq, size_t seq_len, const uint8_t* qual, Sub* sub,
                std::vector<std::vector<double>>& error_matrix, int* tvec, unsigned int* qind);

        // Fill error matrix with values read from the given file.
        // The matrix contains a row for each pair of nucleotides and a column for each possible quality score.
        // Each line of the file corresponds to one row of the matrix, with columns separated by a single space.
        // If the file does not specify enough columns, the last one is repeated.
        static bool fill_error_matrix(std::vector<std::vector<double>>& error_matrix, const std::string& file, short num_scores);

    };


}

#endif //GEFAST_UTILITY_HPP
