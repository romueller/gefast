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

#ifndef GEFAST_PREPROCESSORS_HPP
#define GEFAST_PREPROCESSORS_HPP

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

#include "Base.hpp"
#include "Factories.hpp"
#include "Utility.hpp"

namespace GeFaST {

    /*
     * Reads amplicons from one or more FASTA files into an AmpliconStorage object.
     *
     * All input files are read twice and the sequences are filtered based on the configuration.
     */
    class FastaPreprocessor : public Preprocessor {

        /*
         * Scan the given input file for sequences that satisfy the preprocessor's filters.
         * The sequences passing the filters are accounted for in the DataStatistics object
         * (in order to initialise the AmpliconStorage object appropriately later on).
         */
        void analyse_input_file(const std::string& file, DataStatistics<>& ds, const Configuration& config);

        /*
         * Reread the analysed input and inserts the amplicons satisfying the filters into AmpliconStorage object.
         * "Normalises" the sequences to upper-case letters.
         */
        void read_input_file(const std::string& file, AmpliconStorage& amplicon_storage, const Configuration& config);

    public:
        /*
         * Perform the overall preprocessing step.
         *
         * First, all input files are analysed to determine statistics on the amplicons
         * to be inserted in the AmpliconStorage object.
         * Second, the AmpliconStorage object is initialised based on these statistics.
         * Third, all input files are read a second time to actually insert the amplicons.
         * Finally, amplicons within each pool are sorted by abundance (using the lexicographical order
         * of the headers as the tie-breaker).
         */
        AmpliconStorage* preprocess_inputs(const std::vector<std::string>& files, const Configuration& config) override;

        FastaPreprocessor* clone() const override; // deep-copy clone method

    };


    /*
     * Abstract class for representation of a (dereplicated) FASTQ read.
     */
    struct ReadRepresentation {

        std::string id; // identifier of (first occurrence of) the read
        numSeqs_t abundance; // abundance of the read (over all occurrences)

        /*
         * Initialise read representation with its first occurrence.
         */
        virtual void initialise(const Defline& dl, const std::string& seq, const std::string& quality_scores, const QualityEncoding<>& qe) = 0;

        /*
         * Add abundance and error probabilities (weighted by the abundance) of another occurrence of the read.
         */
        virtual void add_duplicate(const numSeqs_t ab, const std::string& quality_scores, const QualityEncoding<>& qe) = 0;

        /*
         * Get combined error probabilities of dereplicated read.
         */
        virtual std::string get_scores(const QualityEncoding<>& qe) const = 0;

    };

    /*
     * The quality scores of the dereplicated read are set to the maximum quality score (per position)
     * observed for the read.
     */
    struct MaxReadRepresentation : public ReadRepresentation {

        std::string max_scores; // maximum quality score per position (over all occurrences)

        /*
         * Initialise read representation with its first occurrence.
         */
        void initialise(const Defline& dl, const std::string& seq, const std::string& quality_scores, const QualityEncoding<>& qe) override;

        /*
         * Add abundance and error probabilities (weighted by the abundance) of another occurrence of the read.
         */
        void add_duplicate(const numSeqs_t ab, const std::string& quality_scores, const QualityEncoding<>& qe) override;

        /*
         * Get combined error probabilities of dereplicated read.
         */
        std::string get_scores(const QualityEncoding<>& qe) const override;

    };

    /*
     * The quality scores of the dereplicated read correspond to the average error probability (per position)
     * observed for the read.
     */
    struct MeanReadRepresentation : public ReadRepresentation {

        std::vector<numSeqs_t> abundances; // abundance values of all read duplicates
        std::vector<std::string> scores; // quality scores of all read duplicates

        /*
         * Initialise read representation with its first occurrence.
         */
        void initialise(const Defline& dl, const std::string& seq, const std::string& quality_scores, const QualityEncoding<>& qe) override;

        /*
         * Add abundance and error probabilities (weighted by the abundance) of another occurrence of the read.
         */
        void add_duplicate(const numSeqs_t ab, const std::string& quality_scores, const QualityEncoding<>& qe) override;

        /*
         * Get combined error probabilities of dereplicated read.
         */
        std::string get_scores(const QualityEncoding<>& qe) const override;

    };

    /*
     * The quality scores of the dereplicated read are set to the most common quality score (per position)
     * observed for the read.
     */
    struct ModeReadRepresentation : public ReadRepresentation {

        std::vector<std::map<char, numSeqs_t>> score_tallies; // counts of the observed quality scores per position (over all occurrences)

        /*
         * Initialise read representation with its first occurrence.
         */
        void initialise(const Defline& dl, const std::string& seq, const std::string& quality_scores, const QualityEncoding<>& qe) override;

        /*
         * Add abundance and error probabilities (weighted by the abundance) of another occurrence of the read.
         */
        void add_duplicate(const numSeqs_t ab, const std::string& quality_scores, const QualityEncoding<>& qe) override;

        /*
         * Get combined error probabilities of dereplicated read.
         */
        std::string get_scores(const QualityEncoding<>& qe) const override;

    };

    /*
     * The quality scores of the dereplicated read are set to the median quality score (per position)
     * observed for the read.
     */
    struct MedianReadRepresentation : public ReadRepresentation {

        std::vector<std::map<char, numSeqs_t>> score_tallies; // counts of the observed quality scores per position (over all occurrences)

        /*
         * Initialise read representation with its first occurrence.
         */
        void initialise(const Defline& dl, const std::string& seq, const std::string& quality_scores, const QualityEncoding<>& qe) override;

        /*
         * Add abundance and error probabilities (weighted by the abundance) of another occurrence of the read.
         */
        void add_duplicate(const numSeqs_t ab, const std::string& quality_scores, const QualityEncoding<>& qe) override;

        /*
         * Get combined error probabilities of dereplicated read.
         */
        std::string get_scores(const QualityEncoding<>& qe) const override;

    };

    /*
     * The quality scores of the dereplicated read correspond to the product of the error probabilities (per position)
     * observed for the read. In order to stay within range of the used quality encoding, the error probability
     * and, thus, the quality score are capped based on the highest quality score of the encoding.
     */
    struct ProductReadRepresentation : public ReadRepresentation {

        std::vector<float> probability_products; // products of error probability per position (over all occurrences)

        /*
         * Initialise read representation with its first occurrence.
         */
        void initialise(const Defline& dl, const std::string& seq, const std::string& quality_scores, const QualityEncoding<>& qe) override;

        /*
         * Add abundance and error probabilities (weighted by the abundance) of another occurrence of the read.
         */
        void add_duplicate(const numSeqs_t ab, const std::string& quality_scores, const QualityEncoding<>& qe) override;

        /*
         * Get combined error probabilities of dereplicated read.
         */
        std::string get_scores(const QualityEncoding<>& qe) const override;

    };


    /*
     * Reads amplicons from one or more FASTQ files into an AmpliconStorage object.
     *
     * Each input file is read once and the sequences are filtered based on the configuration.
     * The FASTQ entries are dereplicated based on their sequence (to obtain abundance values).
     * During dereplication, the quality scores of sequences grouped together are averaged per position.
     *
     * Whether or not the quality scores are actually stored depends on the used AmpliconStorage class.
     */
    template<typename R>
    class FastqPreprocessor : public Preprocessor {

        /*
         * Scan the given input file for sequences that satisfy the preprocessor's filters.
         * The sequences passing the filters are accounted for in the DataStatistics object
         * (in order to initialise the AmpliconStorage object appropriately later on) and dereplicated.
         */
        void read_input_file(const std::string& file, const Configuration& config, DataStatistics<>& ds,
                std::map<std::string, R>& different_reads) {

            std::ifstream in_stream(file);
            if (!in_stream.good()) {

                std::cerr << " -- WARNING: File '" << file << "' not opened correctly. No sequences are read from it." << std::endl;
                return;

            }

            // set up filters (depending on configuration)
            lenSeqs_t min_length = 0;
            lenSeqs_t max_length = std::numeric_limits<lenSeqs_t>::max();
            switch (config.filter_length) {
                case 1: { // 01 = only max
                    max_length = config.max_length;
                    break;
                }

                case 2: { // 10 = only min
                    min_length = config.min_length;
                    break;
                }

                case 3: { // 11 = min & max
                    min_length = config.min_length;
                    max_length = config.max_length;
                    break;
                }

                default: {
                    // do nothing
                }
            }

            numSeqs_t min_abundance = 0;
            numSeqs_t max_abundance = std::numeric_limits<numSeqs_t>::max();


            Defline dl;
            std::string line, seq;

            while (std::getline(in_stream, line).good()) {

                if (line.empty() || line[0] == ';') continue; // skip empty and comment lines (begin with ';')

                if (line[0] == '@') { // header line

                    dl = Defline::parse_description_line(line, config.separator);

                    std::getline(in_stream, seq); // read sequence
                    upper_case(seq);

                    std::getline(in_stream, line); // read '+' line (ignored)
                    std::getline(in_stream, line); // read quality scores

                    if (check_sequence(seq, config.alphabet, dl.abundance, min_length, max_length, min_abundance, max_abundance)) {

                        auto iter = different_reads.find(seq);
                        if (iter != different_reads.end()) {
                            iter->second.add_duplicate(dl.abundance, line, quality_encoding_);
                        } else {

                            different_reads[seq].initialise(dl, seq, line, quality_encoding_);
                            ds.record(dl.id, seq, dl.abundance);

                        }

                    }

                }

            }

        }

    public:
        explicit FastqPreprocessor(const QualityEncoding<>& qe) : quality_encoding_(qe) {
            // nothing else to do
        }

        explicit FastqPreprocessor(const QualityEncodingOption opt) : quality_encoding_(QualityEncoding<>(opt)) {
            // nothing else to do
        }

        FastqPreprocessor* clone() const override { // deep-copy clone method
            return new FastqPreprocessor(*this);
        }

        /*
         * Perform the overall preprocessing step.
         *
         * First, all input files are analysed and dereplicated to determine statistics on the amplicons
         * to be inserted in the AmpliconStorage object.
         * Second, the AmpliconStorage object is initialised based on these statistics.
         * Third, all dereplicated amplicons are inserted.
         * Finally, amplicons within each pool are sorted by abundance (using the lexicographical order
         * of the headers as the tie-breaker).
         */
        AmpliconStorage* preprocess_inputs(const std::vector<std::string>& files, const Configuration& config) override {

            AmpliconStorage* amplicon_storage = nullptr;
            std::map<std::string, R> different_reads;

            {
                DataStatistics<> ds;

                std::cout << " -- Scanning input files..." << std::endl;
                for (auto& f : files) {
                    read_input_file(f, config, ds, different_reads);
                }

                if (config.filter_abundance > 0) {
                    filter_by_abundance(different_reads, ds, config);
                }

                amplicon_storage = config.build_amplicon_storage(ds);

                std::cout << " -- Scan completed." << std::endl;

            }

            std::cout << " -- Inserting amplicons..." << std::endl;
            for (auto iter = different_reads.begin(); iter != different_reads.end(); iter++) {

                auto& rr = iter->second;
                ExtraInfoQuality ei;
                ei.quality_scores = rr.get_scores(quality_encoding_);
                Defline dl(rr.id, "", rr.abundance);
                amplicon_storage->add_amplicon(dl, iter->first, ei);

            }

            amplicon_storage->finalise();

            std::cout << " -- Construction of amplicon storage completed." << std::endl;

            return amplicon_storage;

        }

    protected:
        FastqPreprocessor() = default;

        /*
         * Remove reads with an abundance that is too high or too low.
         * The filter is set up according to the configuration
         * with no restriction, only a lower bound, only an upper bound, or both.
         */
        void filter_by_abundance(std::map<std::string, R>& different_reads, DataStatistics<>& ds,
                const Configuration& config) {

            numSeqs_t min_abundance = 0;
            numSeqs_t max_abundance = std::numeric_limits<numSeqs_t>::max();
            switch (config.filter_abundance) {
                case 1: { // 01 = only max
                    max_abundance = config.max_abundance;
                    break;
                }

                case 2: { // 10 = only min
                    min_abundance = config.min_abundance;
                    break;
                }

                case 3: { // 11 = min & max
                    min_abundance = config.min_abundance;
                    max_abundance = config.max_abundance;
                    break;
                }

                default: {
                    // do nothing
                }
            }

            for (auto iter = different_reads.begin(); iter != different_reads.end(); ) {

                if (iter->second.abundance < min_abundance || iter->second.abundance > max_abundance) {

                    ds.unrecord(iter->second.id, iter->first, iter->second.abundance);
                    iter = different_reads.erase(iter);

                } else {
                    iter++;
                }

            }

        }

        QualityEncoding<> quality_encoding_; // mapping between used quality scores and error probabilities

    };

}

#endif //GEFAST_PREPROCESSORS_HPP
