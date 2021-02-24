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

#include <iomanip>

#include "../include/AmpliconCollections.hpp"

namespace GeFaST {

    /* === SimpleAmpliconCollection === */

    SimpleAmpliconCollection* SimpleAmpliconCollection::clone() const {
        return new SimpleAmpliconCollection(*this);
    }

    void SimpleAmpliconCollection::emplace_back(const Defline& dl, const std::string& seq, const ExtraInfo& extra) {

        ids_.push_back(dl.id);
        seqs_.push_back(seq);
        abunds_.push_back(dl.abundance);

    }

    numSeqs_t SimpleAmpliconCollection::size() const {
        return ids_.size();
    }

    const char* SimpleAmpliconCollection::id(const numSeqs_t i) const {
        return ids_[i].c_str();
    }

    const char* SimpleAmpliconCollection::seq(const numSeqs_t i) const {
        return seqs_[i].c_str();
    }

    lenSeqs_t SimpleAmpliconCollection::len(const numSeqs_t i) const {
        return seqs_[i].length();
    }

    numSeqs_t SimpleAmpliconCollection::ab(const numSeqs_t i) const {
        return abunds_[i];
    }

    const char* SimpleAmpliconCollection::quals(const numSeqs_t i) const {
        return nullptr;
    }

    unsigned long SimpleAmpliconCollection::qgram_diff(const numSeqs_t i, const numSeqs_t j) const {
        return 0;
    }

    void SimpleAmpliconCollection::sort() {

        std::vector<numSeqs_t> permutation(size());
        std::iota(permutation.begin(), permutation.end(), 0);

        std::sort(permutation.begin(), permutation.end(),
                  [&](const numSeqs_t ampl_a, const numSeqs_t ampl_b) {
                      return (abunds_[ampl_a] > abunds_[ampl_b]) ||
                             ((abunds_[ampl_a] == abunds_[ampl_b]) &&
                              (strcmp(ids_[ampl_a].c_str(), ids_[ampl_b].c_str()) < 0));
                  }
        );

        std::vector<bool> done(size());
        for (numSeqs_t i = 0; i < size(); ++i) {

            if (!done[i]) {

                done[i] = true;
                numSeqs_t prev_j = i;
                numSeqs_t j = permutation[i];

                while (i != j) {

                    std::swap(ids_[prev_j], ids_[j]);
                    std::swap(seqs_[prev_j], seqs_[j]);
                    std::swap(abunds_[prev_j], abunds_[j]);
                    done[j] = true;
                    prev_j = j;
                    j = permutation[j];

                }
            }

        }

    }

    void SimpleAmpliconCollection::print(const std::string& file, const Configuration& config) const {

        std::ofstream out_stream(file, std::ios_base::app);
        std::stringstream str_stream;

        for (numSeqs_t i = 0; i < ids_.size(); i++) {

            str_stream << ">" << ids_[i] << config.separator << abunds_[i] << std::endl << seqs_[i] << std::endl;
            out_stream << str_stream.rdbuf();
            str_stream.str(std::string());

        }

        out_stream.close();

    }


    /* === ArrayAmpliconCollection === */

    ArrayAmpliconCollection::ArrayAmpliconCollection(const numSeqs_t capacity, const unsigned long long capacity_headers,
            const unsigned long long capacity_sequences) {

        capacity_headers_ = capacity_headers + 1;
        headers_ = new char[capacity_headers_];
        next_header_ = headers_;

        capacity_sequences_ = capacity_sequences + 1;
        sequences_ = new char[capacity_sequences_];
        next_sequence_ = sequences_;

        capacity_ = capacity + 1; // +1 to have space for a sentinel after the currently last entry
        header_pointers_ = new char*[capacity_];
        seq_pointers_ = new char*[capacity_];
        abundances_ = new numSeqs_t[capacity_];
        lengths_ = new lenSeqs_t[capacity_];
        size_ = 0;

    }

    ArrayAmpliconCollection::~ArrayAmpliconCollection() {

        delete[] lengths_;
        delete[] abundances_;
        delete[] seq_pointers_;
        delete[] header_pointers_;
        delete[] sequences_;
        delete[] headers_;

    }

    ArrayAmpliconCollection::ArrayAmpliconCollection(const ArrayAmpliconCollection& other) { // copy constructor

        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_headers_ = other.capacity_headers_;
        headers_ = new char[capacity_headers_];
        memcpy(headers_, other.headers_, sizeof(char) * capacity_headers_);
        header_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            header_pointers_[i] = headers_ + (other.header_pointers_[i] - other.headers_);
        }
        next_header_ = headers_ + (other.next_header_ - other.headers_);

        capacity_sequences_ = other.capacity_sequences_;
        sequences_ = new char[capacity_sequences_];
        memcpy(sequences_, other.sequences_, sizeof(char) * capacity_sequences_);
        seq_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            seq_pointers_[i] = sequences_ + (other.seq_pointers_[i] - other.sequences_);
        }
        next_sequence_ = sequences_ + (other.next_sequence_ - other.sequences_);

        abundances_ = new numSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            abundances_[i] = other.abundances_[i];
        }

        lengths_ = new lenSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            lengths_[i] = other.lengths_[i];
        }

    }

    ArrayAmpliconCollection::ArrayAmpliconCollection(ArrayAmpliconCollection&& other) noexcept { // move constructor

        capacity_ = other.capacity_;
        capacity_headers_ = other.capacity_headers_;
        capacity_sequences_ = other.capacity_sequences_;
        size_ = other.size_;
        next_header_ = other.next_header_;
        next_sequence_ = other.next_sequence_;

        headers_ = other.headers_; other.headers_ = nullptr;
        sequences_ = other.sequences_; other.sequences_ = nullptr;
        header_pointers_ = other.header_pointers_; other.header_pointers_ = nullptr;
        seq_pointers_ = other.seq_pointers_; other.seq_pointers_ = nullptr;
        abundances_ = other.abundances_; other.abundances_ = nullptr;
        lengths_ = other.lengths_; other.lengths_ = nullptr;

    }

    ArrayAmpliconCollection& ArrayAmpliconCollection::operator=(const ArrayAmpliconCollection& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] lengths_;
        delete[] abundances_;
        delete[] seq_pointers_;
        delete[] header_pointers_;
        delete[] sequences_;
        delete[] headers_;

        // copy new resources
        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_headers_ = other.capacity_headers_;
        headers_ = new char[capacity_headers_];
        memcpy(headers_, other.headers_, sizeof(char) * capacity_headers_);
        header_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            header_pointers_[i] = headers_ + (other.header_pointers_[i] - other.headers_);
        }
        next_header_ = headers_ + (other.next_header_ - other.headers_);

        capacity_sequences_ = other.capacity_sequences_;
        sequences_ = new char[capacity_sequences_];
        memcpy(sequences_, other.sequences_, sizeof(char) * capacity_sequences_);
        seq_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            seq_pointers_[i] = sequences_ + (other.seq_pointers_[i] - other.sequences_);
        }
        next_sequence_ = sequences_ + (other.next_sequence_ - other.sequences_);

        abundances_ = new numSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            abundances_[i] = other.abundances_[i];
        }

        lengths_ = new lenSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            lengths_[i] = other.lengths_[i];
        }

        return *this;

    }

    ArrayAmpliconCollection& ArrayAmpliconCollection::operator=(ArrayAmpliconCollection&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] lengths_;
        delete[] abundances_;
        delete[] seq_pointers_;
        delete[] header_pointers_;
        delete[] sequences_;
        delete[] headers_;

        // copy / transfer new resources
        capacity_ = other.capacity_;
        capacity_headers_ = other.capacity_headers_;
        capacity_sequences_ = other.capacity_sequences_;
        size_ = other.size_;
        next_header_ = other.next_header_;
        next_sequence_ = other.next_sequence_;

        headers_ = other.headers_; other.headers_ = nullptr;
        sequences_ = other.sequences_; other.sequences_ = nullptr;
        header_pointers_ = other.header_pointers_; other.header_pointers_ = nullptr;
        seq_pointers_ = other.seq_pointers_; other.seq_pointers_ = nullptr;
        abundances_ = other.abundances_; other.abundances_ = nullptr;
        lengths_ = other.lengths_; other.lengths_ = nullptr;

        return *this;

    }

    ArrayAmpliconCollection* ArrayAmpliconCollection::clone() const {
        return new ArrayAmpliconCollection(*this);
    }

    void ArrayAmpliconCollection::emplace_back(const Defline& dl, const std::string& seq, const ExtraInfo& extra) {

        extend_arrays(dl, seq);

        header_pointers_[size_] = next_header_;
        strcpy(next_header_, dl.id.c_str());
        next_header_ += dl.id.length() + 1;

        seq_pointers_[size_] = next_sequence_;
        strcpy(next_sequence_, seq.c_str());
        next_sequence_ += seq.length() + 1;

        abundances_[size_] = dl.abundance;
        lengths_[size_] = seq.length();

        size_++;
        header_pointers_[size_] = next_header_; // sentinel for length calculations
        seq_pointers_[size_] = next_sequence_; // sentinel for length calculations

    }

    numSeqs_t ArrayAmpliconCollection::size() const {
        return size_;
    }

    const char* ArrayAmpliconCollection::id(const numSeqs_t i) const {
        return header_pointers_[i];
    }

    const char* ArrayAmpliconCollection::seq(const numSeqs_t i) const {
        return seq_pointers_[i];
    }

    lenSeqs_t ArrayAmpliconCollection::len(const numSeqs_t i) const {
        return lengths_[i];
    }

    numSeqs_t ArrayAmpliconCollection::ab(const numSeqs_t i) const {
        return abundances_[i];
    }

    const char* ArrayAmpliconCollection::quals(const numSeqs_t i) const {
        return nullptr;
    }

    unsigned long ArrayAmpliconCollection::qgram_diff(const numSeqs_t i, const numSeqs_t j) const {
        return 0;
    }

    void ArrayAmpliconCollection::sort() {

        std::vector<numSeqs_t> permutation(size());
        std::iota(permutation.begin(), permutation.end(), 0);

        std::sort(permutation.begin(), permutation.end(),
                  [&](const numSeqs_t ampl_a, const numSeqs_t ampl_b) {
                      return (abundances_[ampl_a] > abundances_[ampl_b]) ||
                             ((abundances_[ampl_a] == abundances_[ampl_b]) &&
                              (strcmp(header_pointers_[ampl_a], header_pointers_[ampl_b]) < 0));
                  }
        );

        std::vector<bool> done(size());
        for (numSeqs_t i = 0; i < size(); ++i) {

            if (!done[i]) {

                done[i] = true;
                numSeqs_t prev_j = i;
                numSeqs_t j = permutation[i];

                while (i != j) {

                    std::swap(header_pointers_[prev_j], header_pointers_[j]);
                    std::swap(seq_pointers_[prev_j], seq_pointers_[j]);
                    std::swap(abundances_[prev_j], abundances_[j]);
                    std::swap(lengths_[prev_j], lengths_[j]);
                    done[j] = true;
                    prev_j = j;
                    j = permutation[j];

                }
            }

        }

    }

    void ArrayAmpliconCollection::extend_arrays(const Defline& dl, const std::string& seq) {

        unsigned long long len = next_header_ - headers_;
        if (len + dl.id.length() + 1 > capacity_headers_) {

            char* tmp1 = new char[2 * capacity_headers_];
            memcpy(tmp1, headers_, len * sizeof(char));
            char** tmp2 = new char*[capacity_];
            for (auto i = 0; i < size_; i++) {
                tmp2[i] = tmp1 + (header_pointers_[i] - headers_);
            }
            delete[] headers_;
            delete[] header_pointers_;
            headers_ = tmp1;
            next_header_ = headers_ + len;
            header_pointers_ = tmp2;

            capacity_headers_ = std::max(len + dl.id.length() + 1, 2 * capacity_headers_);

        }

        len = next_sequence_ - sequences_;
        if (len + seq.length() + 1 > capacity_sequences_) {

            char* tmp1 = new char[2 * capacity_sequences_];
            memcpy(tmp1, sequences_, len * sizeof(char));
            char** tmp2 = new char*[capacity_];
            for (auto i = 0; i < size_; i++) {
                tmp2[i] = tmp1 + (seq_pointers_[i] - sequences_);
            }
            delete[] sequences_;
            delete[] seq_pointers_;
            sequences_ = tmp1;
            next_sequence_ = sequences_ + len;
            seq_pointers_ = tmp2;

            capacity_sequences_ = std::max(len + seq.length() + 1, 2 * capacity_sequences_);

        }

        if (size_ + 1 == capacity_) {

            char** tmp1 = new char*[2 * capacity_];
            memcpy(tmp1, header_pointers_, capacity_ * sizeof(char*));
            delete[] header_pointers_;
            header_pointers_ = tmp1;

            tmp1 = new char*[2 * capacity_];
            memcpy(tmp1, seq_pointers_, capacity_ * sizeof(char*));
            delete[] seq_pointers_;
            seq_pointers_ = tmp1;

            numSeqs_t* tmp2 = new numSeqs_t[2 * capacity_];
            memcpy(tmp2, abundances_, capacity_ * sizeof(numSeqs_t));
            delete[] abundances_;
            abundances_ = tmp2;

            lenSeqs_t* tmp3 = new lenSeqs_t[2 * capacity_];
            memcpy(tmp3, lengths_, capacity_ * sizeof(lenSeqs_t));
            delete[] lengths_;
            lengths_ = tmp3;

            capacity_ *= 2;

        }

    }

    void ArrayAmpliconCollection::print(const std::string& file, const Configuration& config) const {

        std::ofstream out_stream(file, std::ios_base::app);
        std::stringstream str_stream;

        for (numSeqs_t i = 0; i < size_; i++) {

            str_stream << ">" << header_pointers_[i] << config.separator << abundances_[i] << std::endl
                << seq_pointers_[i] << std::endl;
            out_stream << str_stream.rdbuf();
            str_stream.str(std::string());

        }

        out_stream.close();

    }


    /* === ArrayQualityAmpliconCollection === */

    ArrayQualityAmpliconCollection::ArrayQualityAmpliconCollection() {

        headers_ = new char[1];
        next_header_ = headers_;
        capacity_headers_ = 1;

        sequences_ = new char[1];
        capacity_sequences_ = 1;
        next_sequence_ = sequences_;

        qualities_ = new char[1];
        next_qualities_ = qualities_;

        capacity_ = 1; // +1 to have space for a sentinel after the currently last entry
        header_pointers_ = new char*[capacity_];
        seq_pointers_ = new char*[capacity_];
        abundances_ = new numSeqs_t[capacity_];
        lengths_ = new lenSeqs_t[capacity_];
        size_ = 0;

    }

    ArrayQualityAmpliconCollection::ArrayQualityAmpliconCollection(const numSeqs_t capacity,
            const unsigned long long capacity_headers, const unsigned long long capacity_sequences) {

        headers_ = new char[capacity_headers];
        next_header_ = headers_;
        capacity_headers_ = capacity_headers;

        sequences_ = new char[capacity_sequences];
        capacity_sequences_ = capacity_sequences;
        next_sequence_ = sequences_;

        qualities_ = new char[capacity_sequences];
        next_qualities_ = qualities_;

        capacity_ = capacity + 1; // +1 to have space for a sentinel after the currently last entry
        header_pointers_ = new char*[capacity_];
        seq_pointers_ = new char*[capacity_];
        abundances_ = new numSeqs_t[capacity_];
        lengths_ = new lenSeqs_t[capacity_];
        size_ = 0;

    }

    ArrayQualityAmpliconCollection::~ArrayQualityAmpliconCollection() {

        delete[] lengths_;
        delete[] abundances_;
        delete[] seq_pointers_;
        delete[] header_pointers_;
        delete[] qualities_;
        delete[] sequences_;
        delete[] headers_;

    }

    ArrayQualityAmpliconCollection::ArrayQualityAmpliconCollection(const ArrayQualityAmpliconCollection& other) { // copy constructor

        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_headers_ = other.capacity_headers_;
        headers_ = new char[capacity_headers_];
        memcpy(headers_, other.headers_, sizeof(char) * capacity_headers_);
        header_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            header_pointers_[i] = headers_ + (other.header_pointers_[i] - other.headers_);
        }
        next_header_ = headers_ + (other.next_header_ - other.headers_);

        capacity_sequences_ = other.capacity_sequences_;
        sequences_ = new char[capacity_sequences_];
        memcpy(sequences_, other.sequences_, sizeof(char) * capacity_sequences_);
        seq_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            seq_pointers_[i] = sequences_ + (other.seq_pointers_[i] - other.sequences_);
        }
        next_sequence_ = sequences_ + (other.next_sequence_ - other.sequences_);

        qualities_ = new char[capacity_sequences_];
        memcpy(qualities_, other.qualities_, sizeof(char) * capacity_sequences_);
        next_qualities_ = qualities_ + (other.next_qualities_ - other.qualities_);

        abundances_ = new numSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            abundances_[i] = other.abundances_[i];
        }

        lengths_ = new lenSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            lengths_[i] = other.lengths_[i];
        }

    }

    ArrayQualityAmpliconCollection::ArrayQualityAmpliconCollection(ArrayQualityAmpliconCollection&& other) noexcept { // move constructor

        capacity_ = other.capacity_;
        capacity_headers_ = other.capacity_headers_;
        capacity_sequences_ = other.capacity_sequences_;
        size_ = other.size_;
        next_header_ = other.next_header_;
        next_sequence_ = other.next_sequence_;
        next_qualities_ = other.next_qualities_;

        headers_ = other.headers_; other.headers_ = nullptr;
        sequences_ = other.sequences_; other.sequences_ = nullptr;
        qualities_ = other.qualities_; other.qualities_ = nullptr;
        header_pointers_ = other.header_pointers_; other.header_pointers_ = nullptr;
        seq_pointers_ = other.seq_pointers_; other.seq_pointers_ = nullptr;
        abundances_ = other.abundances_; other.abundances_ = nullptr;
        lengths_ = other.lengths_; other.lengths_ = nullptr;

    }

    ArrayQualityAmpliconCollection& ArrayQualityAmpliconCollection::operator=(const ArrayQualityAmpliconCollection& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] lengths_;
        delete[] abundances_;
        delete[] seq_pointers_;
        delete[] header_pointers_;
        delete[] qualities_;
        delete[] sequences_;
        delete[] headers_;

        // copy new resources
        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_headers_ = other.capacity_headers_;
        headers_ = new char[capacity_headers_];
        memcpy(headers_, other.headers_, sizeof(char) * capacity_headers_);
        header_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            header_pointers_[i] = headers_ + (other.header_pointers_[i] - other.headers_);
        }
        next_header_ = headers_ + (other.next_header_ - other.headers_);

        capacity_sequences_ = other.capacity_sequences_;
        sequences_ = new char[capacity_sequences_];
        memcpy(sequences_, other.sequences_, sizeof(char) * capacity_sequences_);
        seq_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            seq_pointers_[i] = sequences_ + (other.seq_pointers_[i] - other.sequences_);
        }
        next_sequence_ = sequences_ + (other.next_sequence_ - other.sequences_);

        qualities_ = new char[capacity_sequences_];
        memcpy(qualities_, other.qualities_, sizeof(char) * capacity_sequences_);
        next_qualities_ = qualities_ + (other.next_qualities_ - other.qualities_);

        abundances_ = new numSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            abundances_[i] = other.abundances_[i];
        }

        lengths_ = new lenSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            lengths_[i] = other.lengths_[i];
        }

        return *this;

    }

    ArrayQualityAmpliconCollection& ArrayQualityAmpliconCollection::operator=(ArrayQualityAmpliconCollection&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] lengths_;
        delete[] abundances_;
        delete[] seq_pointers_;
        delete[] header_pointers_;
        delete[] qualities_;
        delete[] sequences_;
        delete[] headers_;

        // copy / transfer new resources
        capacity_ = other.capacity_;
        capacity_headers_ = other.capacity_headers_;
        capacity_sequences_ = other.capacity_sequences_;
        size_ = other.size_;
        next_header_ = other.next_header_;
        next_sequence_ = other.next_sequence_;
        next_qualities_ = other.next_qualities_;

        headers_ = other.headers_; other.headers_ = nullptr;
        sequences_ = other.sequences_; other.sequences_ = nullptr;
        qualities_ = other.qualities_; other.qualities_ = nullptr;
        header_pointers_ = other.header_pointers_; other.header_pointers_ = nullptr;
        seq_pointers_ = other.seq_pointers_; other.seq_pointers_ = nullptr;
        abundances_ = other.abundances_; other.abundances_ = nullptr;
        lengths_ = other.lengths_; other.lengths_ = nullptr;

        return *this;

    }

    ArrayQualityAmpliconCollection* ArrayQualityAmpliconCollection::clone() const {
        return new ArrayQualityAmpliconCollection(*this);
    }

    void ArrayQualityAmpliconCollection::emplace_back(const Defline& dl, const std::string& seq, const ExtraInfo& extra) {

        extend_arrays(dl, seq);

        header_pointers_[size_] = next_header_;
        strcpy(next_header_, dl.id.c_str());
        next_header_ += dl.id.length() + 1;

        seq_pointers_[size_] = next_sequence_;
        strcpy(next_sequence_, seq.c_str());
        next_sequence_ += seq.length() + 1;

        auto& ref = (const ExtraInfoQuality&)extra;
        strcpy(next_qualities_, ref.quality_scores.c_str());
        next_qualities_ += ref.quality_scores.length() + 1;

        abundances_[size_] = dl.abundance;
        lengths_[size_] = seq.length();

        size_++;
        header_pointers_[size_] = next_header_; // sentinel for length calculations
        seq_pointers_[size_] = next_sequence_; // sentinel for length calculations

    }

    numSeqs_t ArrayQualityAmpliconCollection::size() const {
        return size_;
    }

    const char* ArrayQualityAmpliconCollection::id(const numSeqs_t i) const {
        return header_pointers_[i];
    }

    const char* ArrayQualityAmpliconCollection::seq(const numSeqs_t i) const {
        return seq_pointers_[i];
    }

    lenSeqs_t ArrayQualityAmpliconCollection::len(const numSeqs_t i) const {
        return lengths_[i];
    }

    numSeqs_t ArrayQualityAmpliconCollection::ab(const numSeqs_t i) const {
        return abundances_[i];
    }

    const char* ArrayQualityAmpliconCollection::quals(const numSeqs_t i) const {
        return qualities_ + (seq_pointers_[i] - sequences_);
    }


    void ArrayQualityAmpliconCollection::sort() {

        std::vector<numSeqs_t> permutation(size());
        std::iota(permutation.begin(), permutation.end(), 0);

        std::sort(permutation.begin(), permutation.end(),
                  [&](const numSeqs_t ampl_a, const numSeqs_t ampl_b) {
                      return (abundances_[ampl_a] > abundances_[ampl_b]) ||
                             ((abundances_[ampl_a] == abundances_[ampl_b]) &&
                              (strcmp(header_pointers_[ampl_a], header_pointers_[ampl_b]) < 0));
                  }
        );

        std::vector<bool> done(size());
        for (numSeqs_t i = 0; i < size(); ++i) {

            if (!done[i]) {

                done[i] = true;
                numSeqs_t prev_j = i;
                numSeqs_t j = permutation[i];

                while (i != j) {

                    std::swap(header_pointers_[prev_j], header_pointers_[j]);
                    std::swap(seq_pointers_[prev_j], seq_pointers_[j]);
                    std::swap(abundances_[prev_j], abundances_[j]);
                    std::swap(lengths_[prev_j], lengths_[j]);
                    done[j] = true;
                    prev_j = j;
                    j = permutation[j];

                }
            }

        }

    }

    void ArrayQualityAmpliconCollection::extend_arrays(const Defline& dl, const std::string& seq) {

        unsigned long long len = next_header_ - headers_;
        if (len + dl.id.length() + 1 > capacity_headers_) {

            char* tmp1 = new char[std::max(len + dl.id.length() + 1, 2 * capacity_headers_)];
            memcpy(tmp1, headers_, len * sizeof(char));
            char** tmp2 = new char*[capacity_];
            for (auto i = 0; i < size_; i++) {
                tmp2[i] = tmp1 + (header_pointers_[i] - headers_);
            }
            delete[] headers_;
            delete[] header_pointers_;
            headers_ = tmp1;
            next_header_ = headers_ + len;
            header_pointers_ = tmp2;

            capacity_headers_ = std::max(len + dl.id.length() + 1, 2 * capacity_headers_);

        }

        len = next_sequence_ - sequences_;
        if (len + seq.length() + 1 > capacity_sequences_) {

            char* tmp1 = new char[std::max(len + seq.length() + 1, 2 * capacity_sequences_)];
            memcpy(tmp1, sequences_, len * sizeof(char));
            char** tmp2 = new char*[capacity_];
            for (auto i = 0; i < size_; i++) {
                tmp2[i] = tmp1 + (seq_pointers_[i] - sequences_);
            }
            delete[] sequences_;
            delete[] seq_pointers_;
            sequences_ = tmp1;
            next_sequence_ = sequences_ + len;
            seq_pointers_ = tmp2;

            tmp1 = new char[std::max(len + seq.length() + 1, 2 * capacity_sequences_)];
            memcpy(tmp1, qualities_, len * sizeof(char));
            delete[] qualities_;
            qualities_ = tmp1;

            capacity_sequences_ = std::max(len + seq.length() + 1, 2 * capacity_sequences_);

        }

        if (size_ + 1 == capacity_) {

            char** tmp1 = new char*[2 * capacity_];
            memcpy(tmp1, header_pointers_, capacity_ * sizeof(char*));
            delete[] header_pointers_;
            header_pointers_ = tmp1;

            tmp1 = new char*[2 * capacity_];
            memcpy(tmp1, seq_pointers_, capacity_ * sizeof(char*));
            delete[] seq_pointers_;
            seq_pointers_ = tmp1;

            numSeqs_t* tmp2 = new numSeqs_t[2 * capacity_];
            memcpy(tmp2, abundances_, capacity_ * sizeof(numSeqs_t));
            delete[] abundances_;
            abundances_ = tmp2;

            lenSeqs_t* tmp3 = new lenSeqs_t[2 * capacity_];
            memcpy(tmp3, lengths_, capacity_ * sizeof(lenSeqs_t));
            delete[] lengths_;
            lengths_ = tmp3;

            capacity_ *= 2;

        }

    }

    void ArrayQualityAmpliconCollection::print(const std::string& file, const Configuration& config) const {

        std::ofstream out_stream(file, std::ios_base::app);
        std::stringstream str_stream;

        for (numSeqs_t i = 0; i < size_; i++) {

            str_stream << "@" << header_pointers_[i] << config.separator << abundances_[i] << std::endl
                << seq_pointers_[i] << std::endl << "+" << std::endl << quals(i) << std::endl;
            out_stream << str_stream.rdbuf();
            str_stream.str(std::string());

        }

        out_stream.close();

    }


    /* === ArrayQgramAmpliconCollection === */

    ArrayQgramAmpliconCollection::ArrayQgramAmpliconCollection(const numSeqs_t capacity,
            const unsigned long long capacity_headers, const unsigned long long capacity_sequences) {

        cpu_features_detect();

        capacity_headers_ = capacity_headers + 1;
        headers_ = new char[capacity_headers_];
        next_header_ = headers_;

        capacity_sequences_ = capacity_sequences + 1;
        sequences_ = new char[capacity_sequences_];
        next_sequence_ = sequences_;

        capacity_ = capacity + 1; // +1 to have space for a sentinel after the currently last entry
        qgrams_ = new unsigned char[capacity_ * QGRAM_VECTOR_BYTES];
        qgram_offsets_ = new numSeqs_t[capacity_];
        header_pointers_ = new char*[capacity_];
        seq_pointers_ = new char*[capacity_];
        abundances_ = new numSeqs_t[capacity_];
        lengths_ = new lenSeqs_t[capacity_];
        size_ = 0;

    }

    ArrayQgramAmpliconCollection::~ArrayQgramAmpliconCollection() {

        delete[] lengths_;
        delete[] abundances_;
        delete[] seq_pointers_;
        delete[] header_pointers_;
        delete[] qgram_offsets_;
        delete[] qgrams_;
        delete[] sequences_;
        delete[] headers_;

    }

    ArrayQgramAmpliconCollection::ArrayQgramAmpliconCollection(const ArrayQgramAmpliconCollection& other) { // copy constructor

        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_headers_ = other.capacity_headers_;
        headers_ = new char[capacity_headers_];
        memcpy(headers_, other.headers_, sizeof(char) * capacity_headers_);
        header_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            header_pointers_[i] = headers_ + (other.header_pointers_[i] - other.headers_);
        }
        next_header_ = headers_ + (other.next_header_ - other.headers_);

        capacity_sequences_ = other.capacity_sequences_;
        sequences_ = new char[capacity_sequences_];
        memcpy(sequences_, other.sequences_, sizeof(char) * capacity_sequences_);
        seq_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            seq_pointers_[i] = sequences_ + (other.seq_pointers_[i] - other.sequences_);
        }
        next_sequence_ = sequences_ + (other.next_sequence_ - other.sequences_);

        qgrams_ = new unsigned char[capacity_ * QGRAM_VECTOR_BYTES];
        memcpy(qgrams_, other.qgrams_, sizeof(unsigned char) * capacity_ * QGRAM_VECTOR_BYTES);
        qgram_offsets_ = new numSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            qgram_offsets_[i] = other.qgram_offsets_[i];
        }

        abundances_ = new numSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            abundances_[i] = other.abundances_[i];
        }

        lengths_ = new lenSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            lengths_[i] = other.lengths_[i];
        }

    }

    ArrayQgramAmpliconCollection::ArrayQgramAmpliconCollection(ArrayQgramAmpliconCollection&& other) noexcept { // move constructor

        capacity_ = other.capacity_;
        capacity_headers_ = other.capacity_headers_;
        capacity_sequences_ = other.capacity_sequences_;
        size_ = other.size_;
        next_header_ = other.next_header_;
        next_sequence_ = other.next_sequence_;

        headers_ = other.headers_; other.headers_ = nullptr;
        sequences_ = other.sequences_; other.sequences_ = nullptr;
        header_pointers_ = other.header_pointers_; other.header_pointers_ = nullptr;
        seq_pointers_ = other.seq_pointers_; other.seq_pointers_ = nullptr;
        qgrams_ = other.qgrams_; other.qgrams_ = nullptr;
        qgram_offsets_ = other.qgram_offsets_; other.qgram_offsets_ = nullptr;
        abundances_ = other.abundances_; other.abundances_ = nullptr;
        lengths_ = other.lengths_; other.lengths_ = nullptr;

    }

    ArrayQgramAmpliconCollection& ArrayQgramAmpliconCollection::operator=(const ArrayQgramAmpliconCollection& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] lengths_;
        delete[] abundances_;
        delete[] seq_pointers_;
        delete[] header_pointers_;
        delete[] qgram_offsets_;
        delete[] qgrams_;
        delete[] sequences_;
        delete[] headers_;

        // copy new resources
        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_headers_ = other.capacity_headers_;
        headers_ = new char[capacity_headers_];
        memcpy(headers_, other.headers_, sizeof(char) * capacity_headers_);
        header_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            header_pointers_[i] = headers_ + (other.header_pointers_[i] - other.headers_);
        }
        next_header_ = headers_ + (other.next_header_ - other.headers_);

        capacity_sequences_ = other.capacity_sequences_;
        sequences_ = new char[capacity_sequences_];
        memcpy(sequences_, other.sequences_, sizeof(char) * capacity_sequences_);
        seq_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            seq_pointers_[i] = sequences_ + (other.seq_pointers_[i] - other.sequences_);
        }
        next_sequence_ = sequences_ + (other.next_sequence_ - other.sequences_);

        qgrams_ = new unsigned char[capacity_ * QGRAM_VECTOR_BYTES];
        memcpy(qgrams_, other.qgrams_, sizeof(unsigned char) * capacity_ * QGRAM_VECTOR_BYTES);
        qgram_offsets_ = new numSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            qgram_offsets_[i] = other.qgram_offsets_[i];
        }

        abundances_ = new numSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            abundances_[i] = other.abundances_[i];
        }

        lengths_ = new lenSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            lengths_[i] = other.lengths_[i];
        }

        return *this;

    }

    ArrayQgramAmpliconCollection& ArrayQgramAmpliconCollection::operator=(ArrayQgramAmpliconCollection&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] lengths_;
        delete[] abundances_;
        delete[] seq_pointers_;
        delete[] header_pointers_;
        delete[] qgram_offsets_;
        delete[] qgrams_;
        delete[] sequences_;
        delete[] headers_;

        // copy / transfer new resources
        capacity_ = other.capacity_;
        capacity_headers_ = other.capacity_headers_;
        capacity_sequences_ = other.capacity_sequences_;
        size_ = other.size_;
        next_header_ = other.next_header_;
        next_sequence_ = other.next_sequence_;

        headers_ = other.headers_; other.headers_ = nullptr;
        sequences_ = other.sequences_; other.sequences_ = nullptr;
        header_pointers_ = other.header_pointers_; other.header_pointers_ = nullptr;
        seq_pointers_ = other.seq_pointers_; other.seq_pointers_ = nullptr;
        qgrams_ = other.qgrams_; other.qgrams_ = nullptr;
        qgram_offsets_ = other.qgram_offsets_; other.qgram_offsets_ = nullptr;
        abundances_ = other.abundances_; other.abundances_ = nullptr;
        lengths_ = other.lengths_; other.lengths_ = nullptr;

        return *this;

    }

    ArrayQgramAmpliconCollection* ArrayQgramAmpliconCollection::clone() const {
        return new ArrayQgramAmpliconCollection(*this);
    }

    void ArrayQgramAmpliconCollection::emplace_back(const Defline& dl, const std::string& seq, const ExtraInfo& extra) {

        extend_arrays(dl, seq);

        header_pointers_[size_] = next_header_;
        strcpy(next_header_, dl.id.c_str());
        next_header_ += dl.id.length() + 1;

        seq_pointers_[size_] = next_sequence_;
        strcpy(next_sequence_, seq.c_str());
        next_sequence_ += seq.length() + 1;

        qgram_offsets_[size_] = size_ * QGRAM_VECTOR_BYTES;
        memset(qgrams_ + qgram_offsets_[size_], 0, QGRAM_VECTOR_BYTES);
        unsigned long qgram = 0;
        unsigned long j = 0;
        while((j < QGRAM_LENGTH - 1) && (j < seq.length())) {

            qgram = (qgram << 2) | (acgtu_map[seq[j]] - 1);
            j++;

        }
        while(j < seq.length()) {

            qgram = (qgram << 2) | (acgtu_map[seq[j]] - 1);
            (qgrams_ + qgram_offsets_[size_])[(qgram >> 3) & (QGRAM_VECTOR_BYTES - 1)] ^= (1 << (qgram & 7));
            j++;

        }

        abundances_[size_] = dl.abundance;
        lengths_[size_] = seq.length();

        size_++;
        header_pointers_[size_] = next_header_; // sentinel for length calculations
        seq_pointers_[size_] = next_sequence_; // sentinel for length calculations

    }

    numSeqs_t ArrayQgramAmpliconCollection::size() const {
        return size_;
    }

    const char* ArrayQgramAmpliconCollection::id(const numSeqs_t i) const {
        return header_pointers_[i];
    }

    const char* ArrayQgramAmpliconCollection::seq(const numSeqs_t i) const {
        return seq_pointers_[i];
    }

    lenSeqs_t ArrayQgramAmpliconCollection::len(const numSeqs_t i) const {
        return lengths_[i];
    }

    numSeqs_t ArrayQgramAmpliconCollection::ab(const numSeqs_t i) const {
        return abundances_[i];
    }

    const char* ArrayQgramAmpliconCollection::quals(const numSeqs_t i) const {
        return nullptr;
    }

    void ArrayQgramAmpliconCollection::sort() {

        std::vector<numSeqs_t> permutation(size());
        std::iota(permutation.begin(), permutation.end(), 0);

        std::sort(permutation.begin(), permutation.end(),
                  [&](const numSeqs_t ampl_a, const numSeqs_t ampl_b) {
                      return (abundances_[ampl_a] > abundances_[ampl_b]) ||
                             ((abundances_[ampl_a] == abundances_[ampl_b]) &&
                              (strcmp(header_pointers_[ampl_a], header_pointers_[ampl_b]) < 0));
                  }
        );

        std::vector<bool> done(size());
        for (numSeqs_t i = 0; i < size(); ++i) {

            if (!done[i]) {

                done[i] = true;
                numSeqs_t prev_j = i;
                numSeqs_t j = permutation[i];

                while (i != j) {

                    std::swap(header_pointers_[prev_j], header_pointers_[j]);
                    std::swap(seq_pointers_[prev_j], seq_pointers_[j]);
                    std::swap(qgram_offsets_[prev_j], qgram_offsets_[j]);
                    std::swap(abundances_[prev_j], abundances_[j]);
                    std::swap(lengths_[prev_j], lengths_[j]);
                    done[j] = true;
                    prev_j = j;
                    j = permutation[j];

                }
            }

        }

    }

    void ArrayQgramAmpliconCollection::extend_arrays(const Defline& dl, const std::string& seq) {

        unsigned long long len = next_header_ - headers_;
        if (len + dl.id.length() + 1 > capacity_headers_) {

            char* tmp1 = new char[2 * capacity_headers_];
            memcpy(tmp1, headers_, len * sizeof(char));
            char** tmp2 = new char*[capacity_];
            for (auto i = 0; i < size_; i++) {
                tmp2[i] = tmp1 + (header_pointers_[i] - headers_);
            }
            delete[] headers_;
            delete[] header_pointers_;
            headers_ = tmp1;
            next_header_ = headers_ + len;
            header_pointers_ = tmp2;

            capacity_headers_ = std::max(len + dl.id.length() + 1, 2 * capacity_headers_);

        }

        len = next_sequence_ - sequences_;
        if (len + seq.length() + 1 > capacity_sequences_) {

            char* tmp1 = new char[2 * capacity_sequences_];
            memcpy(tmp1, sequences_, len * sizeof(char));
            char** tmp2 = new char*[capacity_];
            for (auto i = 0; i < size_; i++) {
                tmp2[i] = tmp1 + (seq_pointers_[i] - sequences_);
            }
            delete[] sequences_;
            delete[] seq_pointers_;
            sequences_ = tmp1;
            next_sequence_ = sequences_ + len;
            seq_pointers_ = tmp2;

            capacity_sequences_ = std::max(len + seq.length() + 1, 2 * capacity_sequences_);

        }

        if (size_ + 1 == capacity_) {

            char** tmp1 = new char*[2 * capacity_];
            memcpy(tmp1, header_pointers_, capacity_ * sizeof(char*));
            delete[] header_pointers_;
            header_pointers_ = tmp1;

            tmp1 = new char*[2 * capacity_];
            memcpy(tmp1, seq_pointers_, capacity_ * sizeof(char*));
            delete[] seq_pointers_;
            seq_pointers_ = tmp1;

            numSeqs_t* tmp2 = new numSeqs_t[2 * capacity_];
            memcpy(tmp2, abundances_, capacity_ * sizeof(numSeqs_t));
            delete[] abundances_;
            abundances_ = tmp2;

            lenSeqs_t* tmp3 = new lenSeqs_t[2 * capacity_];
            memcpy(tmp3, lengths_, capacity_ * sizeof(lenSeqs_t));
            delete[] lengths_;
            lengths_ = tmp3;

            unsigned char* tmp4 = new unsigned char[2 * capacity_ * QGRAM_VECTOR_BYTES];
            numSeqs_t* tmp5 = new numSeqs_t[2 * capacity_];
            memcpy(tmp4, qgrams_, capacity_ * QGRAM_VECTOR_BYTES);
            for (auto i = 0; i < size_; i++) {
                tmp5[i] = qgram_offsets_[i];
            }
            delete[] qgrams_;
            delete[] qgram_offsets_;
            qgrams_ = tmp4;
            qgram_offsets_ = tmp5;

            capacity_ *= 2;

        }

    }

    void ArrayQgramAmpliconCollection::print(const std::string& file, const Configuration& config) const {

        std::ofstream out_stream(file, std::ios_base::app);
        std::stringstream str_stream;

        for (numSeqs_t i = 0; i < size_; i++) {

            str_stream << ">" << header_pointers_[i] << config.separator << abundances_[i] << std::endl
                << seq_pointers_[i] << std::endl;
            out_stream << str_stream.rdbuf();
            str_stream.str(std::string());

        }

        out_stream.close();

    }


    unsigned long ArrayQgramAmpliconCollection::qgram_diff(const numSeqs_t i, const numSeqs_t j) const {

        unsigned long diff_qgrams = compare_qgram_vectors(qgrams_ + qgram_offsets_[i], qgrams_ + qgram_offsets_[j]);
        unsigned long min_diff = (diff_qgrams + 2 * QGRAM_LENGTH - 1) / (2 * QGRAM_LENGTH);
        return min_diff;

    }


    /* === ArrayQgramQualityAmpliconCollection === */

    ArrayQgramQualityAmpliconCollection::ArrayQgramQualityAmpliconCollection(const numSeqs_t capacity,
            const unsigned long long capacity_headers, const unsigned long long capacity_sequences) {

        cpu_features_detect();

        capacity_headers_ = capacity_headers + 1;
        headers_ = new char[capacity_headers_];
        next_header_ = headers_;

        capacity_sequences_ = capacity_sequences + 1;
        sequences_ = new char[capacity_sequences_];
        next_sequence_ = sequences_;

        qualities_ = new char[capacity_sequences_];
        next_qualities_ = qualities_;

        capacity_ = capacity + 1; // +1 to have space for a sentinel after the currently last entry
        qgrams_ = new unsigned char[capacity_ * QGRAM_VECTOR_BYTES];
        qgram_offsets_ = new numSeqs_t[capacity_];
        header_pointers_ = new char*[capacity_];
        seq_pointers_ = new char*[capacity_];
        abundances_ = new numSeqs_t[capacity_];
        lengths_ = new lenSeqs_t[capacity_];
        size_ = 0;

    }

    ArrayQgramQualityAmpliconCollection::~ArrayQgramQualityAmpliconCollection() {

        delete[] lengths_;
        delete[] abundances_;
        delete[] seq_pointers_;
        delete[] header_pointers_;
        delete[] qgram_offsets_;
        delete[] qgrams_;
        delete[] qualities_;
        delete[] sequences_;
        delete[] headers_;

    }

    ArrayQgramQualityAmpliconCollection::ArrayQgramQualityAmpliconCollection(const ArrayQgramQualityAmpliconCollection& other) { // copy constructor

        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_headers_ = other.capacity_headers_;
        headers_ = new char[capacity_headers_];
        memcpy(headers_, other.headers_, sizeof(char) * capacity_headers_);
        header_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            header_pointers_[i] = headers_ + (other.header_pointers_[i] - other.headers_);
        }
        next_header_ = headers_ + (other.next_header_ - other.headers_);

        capacity_sequences_ = other.capacity_sequences_;
        sequences_ = new char[capacity_sequences_];
        memcpy(sequences_, other.sequences_, sizeof(char) * capacity_sequences_);
        seq_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            seq_pointers_[i] = sequences_ + (other.seq_pointers_[i] - other.sequences_);
        }
        next_sequence_ = sequences_ + (other.next_sequence_ - other.sequences_);

        qualities_ = new char[capacity_sequences_];
        memcpy(qualities_, other.qualities_, sizeof(char) * capacity_sequences_);
        next_qualities_ = qualities_ + (other.next_qualities_ - other.qualities_);

        qgrams_ = new unsigned char[capacity_ * QGRAM_VECTOR_BYTES];
        memcpy(qgrams_, other.qgrams_, sizeof(unsigned char) * capacity_ * QGRAM_VECTOR_BYTES);
        qgram_offsets_ = new numSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            qgram_offsets_[i] = other.qgram_offsets_[i];
        }

        abundances_ = new numSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            abundances_[i] = other.abundances_[i];
        }

        lengths_ = new lenSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            lengths_[i] = other.lengths_[i];
        }

    }

    ArrayQgramQualityAmpliconCollection::ArrayQgramQualityAmpliconCollection(ArrayQgramQualityAmpliconCollection&& other) noexcept { // move constructor

        capacity_ = other.capacity_;
        capacity_headers_ = other.capacity_headers_;
        capacity_sequences_ = other.capacity_sequences_;
        size_ = other.size_;
        next_header_ = other.next_header_;
        next_sequence_ = other.next_sequence_;
        next_qualities_ = other.next_qualities_;

        headers_ = other.headers_; other.headers_ = nullptr;
        sequences_ = other.sequences_; other.sequences_ = nullptr;
        qualities_ = other.qualities_; other.qualities_ = nullptr;
        header_pointers_ = other.header_pointers_; other.header_pointers_ = nullptr;
        seq_pointers_ = other.seq_pointers_; other.seq_pointers_ = nullptr;
        qgrams_ = other.qgrams_; other.qgrams_ = nullptr;
        qgram_offsets_ = other.qgram_offsets_; other.qgram_offsets_ = nullptr;
        abundances_ = other.abundances_; other.abundances_ = nullptr;
        lengths_ = other.lengths_; other.lengths_ = nullptr;

    }

    ArrayQgramQualityAmpliconCollection& ArrayQgramQualityAmpliconCollection::operator=(const ArrayQgramQualityAmpliconCollection& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] lengths_;
        delete[] abundances_;
        delete[] seq_pointers_;
        delete[] header_pointers_;
        delete[] qgram_offsets_;
        delete[] qgrams_;
        delete[] qualities_;
        delete[] sequences_;
        delete[] headers_;

        // copy new resources
        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_headers_ = other.capacity_headers_;
        headers_ = new char[capacity_headers_];
        memcpy(headers_, other.headers_, sizeof(char) * capacity_headers_);
        header_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            header_pointers_[i] = headers_ + (other.header_pointers_[i] - other.headers_);
        }
        next_header_ = headers_ + (other.next_header_ - other.headers_);

        capacity_sequences_ = other.capacity_sequences_;
        sequences_ = new char[capacity_sequences_];
        memcpy(sequences_, other.sequences_, sizeof(char) * capacity_sequences_);
        seq_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            seq_pointers_[i] = sequences_ + (other.seq_pointers_[i] - other.sequences_);
        }
        next_sequence_ = sequences_ + (other.next_sequence_ - other.sequences_);

        qualities_ = new char[capacity_sequences_];
        memcpy(qualities_, other.qualities_, sizeof(char) * capacity_sequences_);
        next_qualities_ = qualities_ + (other.next_qualities_ - other.qualities_);

        qgrams_ = new unsigned char[capacity_ * QGRAM_VECTOR_BYTES];
        memcpy(qgrams_, other.qgrams_, sizeof(unsigned char) * capacity_ * QGRAM_VECTOR_BYTES);
        qgram_offsets_ = new numSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            qgram_offsets_[i] = other.qgram_offsets_[i];
        }

        abundances_ = new numSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            abundances_[i] = other.abundances_[i];
        }

        lengths_ = new lenSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            lengths_[i] = other.lengths_[i];
        }

        return *this;

    }

    ArrayQgramQualityAmpliconCollection& ArrayQgramQualityAmpliconCollection::operator=(ArrayQgramQualityAmpliconCollection&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] lengths_;
        delete[] abundances_;
        delete[] seq_pointers_;
        delete[] header_pointers_;
        delete[] qualities_;
        delete[] qgram_offsets_;
        delete[] qgrams_;
        delete[] sequences_;
        delete[] headers_;

        // copy / transfer new resources
        capacity_ = other.capacity_;
        capacity_headers_ = other.capacity_headers_;
        capacity_sequences_ = other.capacity_sequences_;
        size_ = other.size_;
        next_header_ = other.next_header_;
        next_sequence_ = other.next_sequence_;
        next_qualities_ = other.next_qualities_;

        headers_ = other.headers_; other.headers_ = nullptr;
        sequences_ = other.sequences_; other.sequences_ = nullptr;
        qualities_ = other.qualities_; other.qualities_ = nullptr;
        header_pointers_ = other.header_pointers_; other.header_pointers_ = nullptr;
        seq_pointers_ = other.seq_pointers_; other.seq_pointers_ = nullptr;
        qgrams_ = other.qgrams_; other.qgrams_ = nullptr;
        qgram_offsets_ = other.qgram_offsets_; other.qgram_offsets_ = nullptr;
        abundances_ = other.abundances_; other.abundances_ = nullptr;
        lengths_ = other.lengths_; other.lengths_ = nullptr;

        return *this;

    }

    ArrayQgramQualityAmpliconCollection* ArrayQgramQualityAmpliconCollection::clone() const {
        return new ArrayQgramQualityAmpliconCollection(*this);
    }

    void ArrayQgramQualityAmpliconCollection::emplace_back(const Defline& dl, const std::string& seq, const ExtraInfo& extra) {

        extend_arrays(dl, seq);

        header_pointers_[size_] = next_header_;
        strcpy(next_header_, dl.id.c_str());
        next_header_ += dl.id.length() + 1;

        seq_pointers_[size_] = next_sequence_;
        strcpy(next_sequence_, seq.c_str());
        next_sequence_ += seq.length() + 1;

        auto& ref = (const ExtraInfoQuality&)extra;
        strcpy(next_qualities_, ref.quality_scores.c_str());
        next_qualities_ += ref.quality_scores.length() + 1;

        qgram_offsets_[size_] = size_ * QGRAM_VECTOR_BYTES;
        memset(qgrams_ + qgram_offsets_[size_], 0, QGRAM_VECTOR_BYTES);
        unsigned long qgram = 0;
        unsigned long j = 0;
        while((j < QGRAM_LENGTH - 1) && (j < seq.length())) {

            qgram = (qgram << 2) | (acgtu_map[seq[j]] - 1);
            j++;

        }
        while(j < seq.length()) {

            qgram = (qgram << 2) | (acgtu_map[seq[j]] - 1);
            (qgrams_ + qgram_offsets_[size_])[(qgram >> 3) & (QGRAM_VECTOR_BYTES - 1)] ^= (1 << (qgram & 7));
            j++;

        }

        abundances_[size_] = dl.abundance;
        lengths_[size_] = seq.length();

        size_++;
        header_pointers_[size_] = next_header_; // sentinel for length calculations
        seq_pointers_[size_] = next_sequence_; // sentinel for length calculations

    }

    numSeqs_t ArrayQgramQualityAmpliconCollection::size() const {
        return size_;
    }

    const char* ArrayQgramQualityAmpliconCollection::id(const numSeqs_t i) const {
        return header_pointers_[i];
    }

    const char* ArrayQgramQualityAmpliconCollection::seq(const numSeqs_t i) const {
        return seq_pointers_[i];
    }

    lenSeqs_t ArrayQgramQualityAmpliconCollection::len(const numSeqs_t i) const {
        return lengths_[i];
    }

    numSeqs_t ArrayQgramQualityAmpliconCollection::ab(const numSeqs_t i) const {
        return abundances_[i];
    }

    const char* ArrayQgramQualityAmpliconCollection::quals(const numSeqs_t i) const {
        return qualities_ + (seq_pointers_[i] - sequences_);
    }


    void ArrayQgramQualityAmpliconCollection::sort() {

        std::vector<numSeqs_t> permutation(size());
        std::iota(permutation.begin(), permutation.end(), 0);

        std::sort(permutation.begin(), permutation.end(),
                  [&](const numSeqs_t ampl_a, const numSeqs_t ampl_b) {
                      return (abundances_[ampl_a] > abundances_[ampl_b]) ||
                             ((abundances_[ampl_a] == abundances_[ampl_b]) &&
                              (strcmp(header_pointers_[ampl_a], header_pointers_[ampl_b]) < 0));
                  }
        );

        std::vector<bool> done(size());
        for (numSeqs_t i = 0; i < size(); ++i) {

            if (!done[i]) {

                done[i] = true;
                numSeqs_t prev_j = i;
                numSeqs_t j = permutation[i];

                while (i != j) {

                    std::swap(header_pointers_[prev_j], header_pointers_[j]);
                    std::swap(seq_pointers_[prev_j], seq_pointers_[j]);
                    std::swap(qgram_offsets_[prev_j], qgram_offsets_[j]);
                    std::swap(abundances_[prev_j], abundances_[j]);
                    std::swap(lengths_[prev_j], lengths_[j]);
                    done[j] = true;
                    prev_j = j;
                    j = permutation[j];

                }
            }

        }

    }

    void ArrayQgramQualityAmpliconCollection::extend_arrays(const Defline& dl, const std::string& seq) {

        unsigned long long len = next_header_ - headers_;
        if (len + dl.id.length() + 1 > capacity_headers_) {

            char* tmp1 = new char[2 * capacity_headers_];
            memcpy(tmp1, headers_, len * sizeof(char));
            char** tmp2 = new char*[capacity_];
            for (auto i = 0; i < size_; i++) {
                tmp2[i] = tmp1 + (header_pointers_[i] - headers_);
            }
            delete[] headers_;
            delete[] header_pointers_;
            headers_ = tmp1;
            next_header_ = headers_ + len;
            header_pointers_ = tmp2;

            capacity_headers_ = std::max(len + dl.id.length() + 1, 2 * capacity_headers_);

        }

        len = next_sequence_ - sequences_;
        if (len + seq.length() + 1 > capacity_sequences_) {

            char* tmp1 = new char[2 * capacity_sequences_];
            memcpy(tmp1, sequences_, len * sizeof(char));
            char** tmp2 = new char*[capacity_];
            for (auto i = 0; i < size_; i++) {
                tmp2[i] = tmp1 + (seq_pointers_[i] - sequences_);
            }
            delete[] sequences_;
            delete[] seq_pointers_;
            sequences_ = tmp1;
            next_sequence_ = sequences_ + len;
            seq_pointers_ = tmp2;

            tmp1 = new char[std::max(len + seq.length() + 1, 2 * capacity_sequences_)];
            memcpy(tmp1, qualities_, len * sizeof(char));
            delete[] qualities_;
            qualities_ = tmp1;

            capacity_sequences_ = std::max(len + seq.length() + 1, 2 * capacity_sequences_);

        }

        if (size_ + 1 == capacity_) {

            char** tmp1 = new char*[2 * capacity_];
            memcpy(tmp1, header_pointers_, capacity_ * sizeof(char*));
            delete[] header_pointers_;
            header_pointers_ = tmp1;

            tmp1 = new char*[2 * capacity_];
            memcpy(tmp1, seq_pointers_, capacity_ * sizeof(char*));
            delete[] seq_pointers_;
            seq_pointers_ = tmp1;

            numSeqs_t* tmp2 = new numSeqs_t[2 * capacity_];
            memcpy(tmp2, abundances_, capacity_ * sizeof(numSeqs_t));
            delete[] abundances_;
            abundances_ = tmp2;

            lenSeqs_t* tmp3 = new lenSeqs_t[2 * capacity_];
            memcpy(tmp3, lengths_, capacity_ * sizeof(lenSeqs_t));
            delete[] lengths_;
            lengths_ = tmp3;

            unsigned char* tmp4 = new unsigned char[2 * capacity_ * QGRAM_VECTOR_BYTES];
            numSeqs_t* tmp5 = new numSeqs_t[2 * capacity_];
            memcpy(tmp4, qgrams_, capacity_ * QGRAM_VECTOR_BYTES);
            for (auto i = 0; i < size_; i++) {
                tmp5[i] = qgram_offsets_[i];
            }
            delete[] qgrams_;
            delete[] qgram_offsets_;
            qgrams_ = tmp4;
            qgram_offsets_ = tmp5;

            capacity_ *= 2;

        }

    }

    void ArrayQgramQualityAmpliconCollection::print(const std::string& file, const Configuration& config) const {

        std::ofstream out_stream(file, std::ios_base::app);
        std::stringstream str_stream;

        for (numSeqs_t i = 0; i < size_; i++) {

            str_stream << "@" << header_pointers_[i] << config.separator << abundances_[i] << std::endl
                << seq_pointers_[i] << std::endl << "+" << std::endl << quals(i) << std::endl;
            out_stream << str_stream.rdbuf();
            str_stream.str(std::string());

        }

        out_stream.close();

    }


    unsigned long ArrayQgramQualityAmpliconCollection::qgram_diff(const numSeqs_t i, const numSeqs_t j) const {

        unsigned long diff_qgrams = compare_qgram_vectors(qgrams_ + qgram_offsets_[i], qgrams_ + qgram_offsets_[j]);
        unsigned long min_diff = (diff_qgrams + 2 * QGRAM_LENGTH - 1) / (2 * QGRAM_LENGTH);
        return min_diff;

    }

}