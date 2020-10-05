/*
 * GeFaST
 *
 * Copyright (C) 2016 - 2020 Robert Mueller
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
     * Reads amplicons from one or more FASTQ files into an AmpliconStorage object.
     *
     * Each input file is read once and the sequences are filtered based on the configuration.
     * The FASTQ entries are dereplicated based on their sequence (to obtain abundance values).
     * During dereplication, the quality scores of sequences grouped together are averaged per position.
     *
     * Whether or not the quality scores are actually stored depends on the used AmpliconStorage class.
     */
    class FastqPreprocessor : public Preprocessor {

        /*
         * Representation of a (dereplicated) FASTQ read.
         */
        struct ReadRepresentation {

            std::string id; // identifier of (first occurrence of) the read
            numSeqs_t abundance; // abundance of the read (over all occurrences)
            std::vector<float> probability_sums; // sum of error probability per position (over all occurrences)

            ReadRepresentation() = default;

            /*
             * Initialise read representation with its first occurrence.
             */
            ReadRepresentation(const Defline& dl, const std::string& seq, const std::vector<float>& quality_scores);

            /*
             * Add abundance and error probabilities (weighted by the abundance) of another occurrence of the read.
             */
            void add_duplicate(const numSeqs_t ab, const std::vector<float>& quality_scores);

        };

        /*
         * Scan the given input file for sequences that satisfy the preprocessor's filters.
         * The sequences passing the filters are accounted for in the DataStatistics object
         * (in order to initialise the AmpliconStorage object appropriately later on) and dereplicated.
         */
        void read_input_file(const std::string& file, const Configuration& config, DataStatistics<>& ds,
                std::map<std::string, ReadRepresentation>& different_reads);

    public:
        explicit FastqPreprocessor(const QualityEncoding<>& qe);

        explicit FastqPreprocessor(const QualityEncodingOption opt);

        FastqPreprocessor* clone() const override; // deep-copy clone method

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
        AmpliconStorage* preprocess_inputs(const std::vector<std::string>& files, const Configuration& config) override;

    protected:
        FastqPreprocessor() = default;

        /*
         * Remove reads with an abundance that is too high or too low.
         * The filter is set up according to the configuration
         * with no restriction, only a lower bound, only an upper bound, or both.
         */
        void filter_by_abundance(std::map<std::string, ReadRepresentation>& different_reads, DataStatistics<>& ds,
                const Configuration& config);

        QualityEncoding<> quality_encoding_; // mapping between used quality scores and error probabilities

    };

}

#endif //GEFAST_PREPROCESSORS_HPP
