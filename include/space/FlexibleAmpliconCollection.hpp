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

#ifndef GEFAST_FLEXIBLEAMPLICONCOLLECTION_HPP
#define GEFAST_FLEXIBLEAMPLICONCOLLECTION_HPP

#include <fstream>
#include <numeric>
#include <sstream>

#include "../Base.hpp"
#include "../modes/SpaceLevenshteinMode.hpp"
#include "Basics.hpp"
#include "dacs/Utility.hpp"

namespace GeFaST {

    /*
     * Flexible amplicon collection storing the headers / identifiers, sequences, lengths, abundances
     * and q-gram profiles in separate, exchangeable components for the investigation of
     * space-efficient data structures.
     *
     * The components are not built during by the constructor of FlexibleAmpliconCollection,
     * but are constructed outside and then added to the collection through setter methods.
     * The FlexibleAmpliconCollection takes ownership over its components.
     */
    class FlexibleAmpliconCollection : public AmpliconCollection {

    public:
        FlexibleAmpliconCollection(const numSeqs_t capacity) : capacity_(capacity), size_(0) {

            cpu_features_detect();

            // components are initialised by setter methods

        }

        ~FlexibleAmpliconCollection() override {

            delete lengths_;
            delete abundances_;
            delete qgrams_;
            delete sequences_;
            delete headers_;

        }

        FlexibleAmpliconCollection(const FlexibleAmpliconCollection& other)
                : capacity_(other.capacity_), size_(other.size_), headers_(other.headers_->clone()), sequences_(other.sequences_->clone()),
                  abundances_(other.abundances_->clone()), lengths_(other.lengths_->clone()), qgrams_(other.qgrams_->clone()) { // copy constructor

            // nothing else to do

        }

        FlexibleAmpliconCollection(FlexibleAmpliconCollection&& other) noexcept { // move constructor

            capacity_ = other.capacity_;
            size_ = other.size_;

            headers_ = other.headers_; other.headers_ = nullptr;
            sequences_ = other.sequences_; other.sequences_ = nullptr;
            qgrams_ = other.qgrams_; other.qgrams_ = nullptr;
            abundances_ = other.abundances_; other.abundances_ = nullptr;
            lengths_ = other.lengths_; other.lengths_ = nullptr;

        }

        FlexibleAmpliconCollection& operator=(const FlexibleAmpliconCollection& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete lengths_;
            delete abundances_;
            delete qgrams_;
            delete sequences_;
            delete headers_;

            // copy new resources
            capacity_ = other.capacity_;
            size_ = other.size_;

            headers_ = other.headers_->clone();
            sequences_ = other.sequences_->clone();
            qgrams_ = other.qgrams_->clone();
            abundances_ = other.abundances_->clone();
            lengths_ = other.lengths_->clone();

            return *this;

        }

        FlexibleAmpliconCollection& operator=(FlexibleAmpliconCollection&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete lengths_;
            delete abundances_;
            delete qgrams_;
            delete sequences_;
            delete headers_;

            // copy / transfer new resources
            capacity_ = other.capacity_;
            size_ = other.size_;

            headers_ = other.headers_; other.headers_ = nullptr;
            sequences_ = other.sequences_; other.sequences_ = nullptr;
            qgrams_ = other.qgrams_; other.qgrams_ = nullptr;
            abundances_ = other.abundances_; other.abundances_ = nullptr;
            lengths_ = other.lengths_; other.lengths_ = nullptr;

            return *this;

        }

        FlexibleAmpliconCollection* clone() const override { // deep-copy clone method
            return new FlexibleAmpliconCollection(*this);
        }


        /* component setters for headers, sequences, lengths, abundances and q-gram profiles */

        void set_header_component(StringComponent* c) {

            delete headers_;
            headers_ = c;

        }

        void set_sequence_component(StringComponent* c) {

            delete sequences_;
            sequences_ = c;

        }

        void set_abundance_component(NumericComponent<numSeqs_t>* c) {

            delete abundances_;
            abundances_ = c;

        }

        void set_length_component(NumericComponent<lenSeqs_t>* c) {

            delete lengths_;
            lengths_ = c;

        }

        void set_qgram_component(QgramComponent* c) {

            delete qgrams_;
            qgrams_ = c;

        }


        /* usual AmpliconCollection interface */

        void emplace_back(const Defline& dl, const std::string& seq, const ExtraInfo& extra) override {

            extend_arrays(dl, seq);

            headers_->push_back(dl.id);
            sequences_->push_back(seq);
            qgrams_->push_back(seq);
            abundances_->push_back(dl.abundance);
            lengths_->push_back(seq.length());

            size_++;

        }

        numSeqs_t size() const override {
            return size_;
        }

        inline std::unique_ptr<SequenceWrapper> id(const numSeqs_t i) const override {
            return headers_->get(i, 0); // second parameter never used (needed due to common interface with sequences)
        }

        inline std::string id_str(const numSeqs_t i) const override {
            return headers_->get_str(i, 0); // second parameter never used (needed due to common interface with sequences)
        }

        inline std::unique_ptr<SequenceWrapper> seq(const numSeqs_t i) const override {
            return sequences_->get(i, lengths_->get(i));
        }

        inline std::string seq_str(const numSeqs_t i) const override {
            return sequences_->get_str(i, lengths_->get(i));
        }

        inline lenSeqs_t len(const numSeqs_t i) const override {
            return lengths_->get(i);
        }

        inline numSeqs_t ab(const numSeqs_t i) const override {
            return abundances_->get(i);
        }

        inline const char* quals(const numSeqs_t i) const override {
            return nullptr;
        }

        /*
         * Compute lower bound on the q-gram distance between two amplicons.
         *
         * Based on the ideas from "Approximate string-matching with q-grams and maximal matches"
         * (Ukkonen, 1992, https://core.ac.uk/download/pdf/82613520.pdf),
         * as implemented in Swarm (https://github.com/torognes/swarm).
         *
         * Only considers whether a q-gram has an even (including 0) or odd number of occurrences,
         * instead of using actual counts of the q-grams.
         */
        unsigned long qgram_diff(const numSeqs_t i, const numSeqs_t j) const override {

            unsigned long diff_qgrams = compare_qgram_vectors(qgrams_->get(i), qgrams_->get(j));
            unsigned long min_diff = (diff_qgrams + 2 * QGRAM_LENGTH - 1) / (2 * QGRAM_LENGTH);
            return min_diff;

        }

        /*
         * Sort the amplicons by abundance (decreasing) and
         * use the lexicographical rank (ascending) of the identifiers as the tie-breaker.
         */
        void sort() {

			headers_->finalise();
			sequences_->finalise();
			abundances_->finalise();
			lengths_->finalise();
			qgrams_->finalise();

            std::vector<numSeqs_t> permutation(size());
            std::iota(permutation.begin(), permutation.end(), 0);

            // second parameter of headers_->get() never used (needed due to common interface with sequences)
            std::sort(permutation.begin(), permutation.end(),
                      [&](const numSeqs_t ampl_a, const numSeqs_t ampl_b) {
                          return (abundances_->get(ampl_a) > abundances_->get(ampl_b)) ||
                                 ((abundances_->get(ampl_a) == abundances_->get(ampl_b)) &&
                                  (*headers_->get(ampl_a, 0) < *headers_->get(ampl_b, 0)));
                      }
            );

			headers_->sort(permutation);
			sequences_->sort(permutation);
			abundances_->sort(permutation);
			lengths_->sort(permutation);
			qgrams_->sort(permutation);

			headers_->finalise();
			sequences_->finalise();
			abundances_->finalise();
			lengths_->finalise();
			qgrams_->finalise();

        }

        /*
         * Append amplicon to the collection, extending its capacity if necessary.
         *
         * Ideally, the amplicon collection is initialised with sufficient space (see the constructor)
         * but if the number of amplicons gets too large or the space for the identifiers resp. sequences
         * becomes insufficient, the capacity of the array is doubled.
         */
        void extend_arrays(const Defline& dl, const std::string& seq) {

            headers_->extend(dl.id);
            sequences_->extend(seq);
            qgrams_->extend();
            abundances_->extend();
            lengths_->extend();

            if (size_ == capacity_) capacity_ *= 2;

        }

        /*
         * Append the amplicons in the stored order in FASTA format to the specified file.
         */
        void print(const std::string& file, const Configuration& config) const {

            std::ofstream out_stream(file, std::ios_base::app);
            std::stringstream str_stream;

            for (numSeqs_t i = 0; i < size_; i++) {

                // second parameter of headers_->get() never used (needed due to common interface with sequences)
                str_stream << ">" << *headers_->get(i, 0) << config.separator << abundances_->get(i) << std::endl
                           << *sequences_->get(i, lengths_->get(i)) << std::endl;
                out_stream << str_stream.rdbuf();
                str_stream.str(std::string());

            }

            out_stream.close();

        }

        size_t size_in_bytes() const override {
            return sizeof(FlexibleAmpliconCollection) // FlexibleAmpliconCollection itself
                + ((headers_ != nullptr) ? headers_->size_in_bytes() : 0) // headers_
                + ((sequences_ != nullptr) ? sequences_->size_in_bytes() : 0) // sequences_
                + ((abundances_ != nullptr) ? abundances_->size_in_bytes() : 0) // abundances_
                + ((lengths_ != nullptr) ? lengths_->size_in_bytes() : 0) // lengths_
                + ((qgrams_ != nullptr) ? qgrams_->size_in_bytes() : 0); // qgrams_
        }

        void show_memory(numSeqs_t pid) const override {

            std::cout << "##################################################" << std::endl;
            std::cout << "# FlexibleAmpliconCollection " << std::endl;
            std::cout << "#  - pid: " << pid << std::endl;
            std::cout << "#  - size: " << size() << std::endl;
            std::cout << "# " << std::endl;
            std::cout << "# sizeof(FlexibleAmpliconCollection): " << sizeof(FlexibleAmpliconCollection) << " bytes" << std::endl;
            std::cout << "# sizeof(StringComponent*): " << sizeof(StringComponent*) << " bytes" << std::endl;
            std::cout << "# sizeof(NumericComponent<numSeqs_t>*): " << sizeof(NumericComponent<numSeqs_t>*) << " bytes" << std::endl;
            std::cout << "# sizeof(NumericComponent<lenSeqs_t>*): " << sizeof(NumericComponent<lenSeqs_t>*) << " bytes" << std::endl;
            std::cout << "# sizeof(QgramComponent*): " << sizeof(QgramComponent*) << " bytes" << std::endl;
            std::cout << "# sizeof(numSeqs_t): " << sizeof(numSeqs_t) << " bytes" << std::endl;
            std::cout << "# " << std::endl;
            std::cout << "# Total size: " << size_in_bytes() << " bytes" << std::endl;
            std::cout << "# headers_: " << ((headers_ != nullptr) ? headers_->size_in_bytes() : 0) << " bytes" << std::endl;
            if (headers_ != nullptr) headers_->show_memory();
            std::cout << "# sequences_: " << ((sequences_ != nullptr) ? sequences_->size_in_bytes() : 0) << " bytes" << std::endl;
            if (sequences_ != nullptr) sequences_->show_memory();
            std::cout << "# abundances_: " << ((abundances_ != nullptr) ? abundances_->size_in_bytes() : 0) << " bytes" << std::endl;
            if (abundances_ != nullptr) abundances_->show_memory();
            std::cout << "# lengths_: " << ((lengths_ != nullptr) ? lengths_->size_in_bytes() : 0) << " bytes" << std::endl;
            if (lengths_ != nullptr) lengths_->show_memory();
            std::cout << "# qgrams_: " << ((qgrams_ != nullptr) ? qgrams_->size_in_bytes() : 0) << " bytes" << std::endl;
            if (qgrams_ != nullptr) qgrams_->show_memory();
            std::cout << "##################################################" << std::endl;

        }

    protected:
        FlexibleAmpliconCollection() = default;

        StringComponent* headers_ = nullptr; // container storing all identifiers
        StringComponent* sequences_ = nullptr; // container storing all sequences
        NumericComponent<numSeqs_t>* abundances_ = nullptr; // container storing all abundances
        NumericComponent<lenSeqs_t>* lengths_ = nullptr; // container storing all lengths
        QgramComponent* qgrams_ = nullptr; // container storing all q-gram profiles

        numSeqs_t size_ = 0; // current number of amplicons in the collection
        numSeqs_t capacity_ = 0; // capacity of the collection (in terms of the number of amplicons)

    };

}


#endif //GEFAST_FLEXIBLEAMPLICONCOLLECTION_HPP