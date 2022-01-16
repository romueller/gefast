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

#include "../../include/space/NumericArrays.hpp"

namespace GeFaST {

    /* === QgramsArray === */

    QgramsArray::QgramsArray(const unsigned long long capacity) {

        size_ = 0;
        capacity_ = capacity;
        qgrams_ = new unsigned char[capacity_ * QGRAM_VECTOR_BYTES];
        qgram_offsets_ = new numSeqs_t[capacity_];

    }

    QgramsArray::~QgramsArray() {

        delete[] qgram_offsets_;
        delete[] qgrams_;

    }

    QgramsArray::QgramsArray(const QgramsArray& other) { // copy constructor

        capacity_ = other.capacity_;
        size_ = other.size_;

        qgrams_ = new unsigned char[capacity_ * QGRAM_VECTOR_BYTES];
        memcpy(qgrams_, other.qgrams_, sizeof(unsigned char) * capacity_ * QGRAM_VECTOR_BYTES);
        qgram_offsets_ = new numSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            qgram_offsets_[i] = other.qgram_offsets_[i];
        }

    }

    QgramsArray::QgramsArray(QgramsArray&& other) noexcept { // move constructor

        capacity_ = other.capacity_;
        size_ = other.size_;

        qgrams_ = other.qgrams_; other.qgrams_ = nullptr;
        qgram_offsets_ = other.qgram_offsets_; other.qgram_offsets_ = nullptr;

    }

    QgramsArray& QgramsArray::operator=(const QgramsArray& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] qgram_offsets_;
        delete[] qgrams_;

        // copy new resources
        capacity_ = other.capacity_;
        size_ = other.size_;

        qgrams_ = new unsigned char[capacity_ * QGRAM_VECTOR_BYTES];
        memcpy(qgrams_, other.qgrams_, sizeof(unsigned char) * capacity_ * QGRAM_VECTOR_BYTES);
        qgram_offsets_ = new numSeqs_t[capacity_];
        for (auto i = 0; i < size_; i++) {
            qgram_offsets_[i] = other.qgram_offsets_[i];
        }

        return *this;

    }

    QgramsArray& QgramsArray::operator=(QgramsArray&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] qgram_offsets_;
        delete[] qgrams_;

        // copy / transfer new resources
        capacity_ = other.capacity_;
        size_ = other.size_;

        qgrams_ = other.qgrams_; other.qgrams_ = nullptr;
        qgram_offsets_ = other.qgram_offsets_; other.qgram_offsets_ = nullptr;

        return *this;

    }

    QgramsArray* QgramsArray::clone() const { // deep-copy clone method
        return new QgramsArray(*this);
    }

    void QgramsArray::push_back(const std::string& seq) {

        extend();

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

        size_++;

    }

    unsigned char* QgramsArray::get(const numSeqs_t i) const {
        return qgrams_ + qgram_offsets_[i];
    }

    void QgramsArray::extend() {

        if (size_ == capacity_) {

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

    void QgramsArray::finalise() {
        // nothing to do
    }

    void QgramsArray::sort(const std::vector<numSeqs_t>& permutation) {

        std::vector<bool> done(size_);
        for (numSeqs_t i = 0; i < size_; ++i) {

            if (!done[i]) {

                done[i] = true;
                numSeqs_t prev_j = i;
                numSeqs_t j = permutation[i];

                while (i != j) {

                    // swap_values(prev_j, j);
                    std::swap(qgram_offsets_[prev_j], qgram_offsets_[j]);

                    done[j] = true;
                    prev_j = j;
                    j = permutation[j];

                }

            }

        }

    }

    size_t QgramsArray::size_in_bytes() const {
        return sizeof(QgramsArray) // QgramsArray itself
            + sizeof(unsigned char) * capacity_ * QGRAM_VECTOR_BYTES // qgrams_
            + sizeof(numSeqs_t) * capacity_; // qgram_offsets_
    }

    void QgramsArray::show_memory() const {

        std::cout << "#  * QgramsArray " << std::endl;
        std::cout << "#  * " << std::endl;
        std::cout << "#  * sizeof(QgramsArray): " << sizeof(QgramsArray) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char): " << sizeof(unsigned char) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char*): " << sizeof(unsigned char*) << " bytes" << std::endl;
        std::cout << "#  * sizeof(numSeqs_t): " << sizeof(numSeqs_t) << " bytes" << std::endl;
        std::cout << "#  * sizeof(numSeqs_t*): " << sizeof(numSeqs_t*) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned long long): " << sizeof(unsigned long long) << " bytes" << std::endl;
        std::cout << "#  * " << std::endl;
        std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
        std::cout << "#  * qgrams_: " << (sizeof(unsigned char) * capacity_ * QGRAM_VECTOR_BYTES) << " bytes" << std::endl;
        std::cout << "#  * qgram_offsets_: " << (sizeof(numSeqs_t) * capacity_) << " bytes" << std::endl;

    }


    /* === QgramsDacsArray === */

    QgramsDacsArray::QgramsDacsArray(const unsigned long long capacity, const ChunkLengths& config) {

        size_ = 0;
        capacity_ = capacity;
        qgrams_ = new unsigned char[capacity_ * QGRAM_VECTOR_BYTES];
        qgram_offsets_ = DacsCapacityArray<dacs_t, builder_t, numSeqs_t>(config.qgram_chunk, capacity, config.qgram_level,
                                                                         config.qgram_extend, config.qgram_opt);

    }

    QgramsDacsArray::~QgramsDacsArray() {
        delete[] qgrams_;
    }

    QgramsDacsArray::QgramsDacsArray(const QgramsDacsArray& other) { // copy constructor

        capacity_ = other.capacity_;
        size_ = other.size_;

        qgrams_ = new unsigned char[capacity_ * QGRAM_VECTOR_BYTES];
        memcpy(qgrams_, other.qgrams_, sizeof(unsigned char) * capacity_ * QGRAM_VECTOR_BYTES);
        qgram_offsets_ = other.qgram_offsets_;

    }

    QgramsDacsArray::QgramsDacsArray(QgramsDacsArray&& other) noexcept { // move constructor

        capacity_ = other.capacity_;
        size_ = other.size_;

        qgrams_ = other.qgrams_; other.qgrams_ = nullptr;
        qgram_offsets_ = std::move(other.qgram_offsets_);

    }

    QgramsDacsArray& QgramsDacsArray::operator=(const QgramsDacsArray& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] qgrams_;

        // copy new resources
        capacity_ = other.capacity_;
        size_ = other.size_;

        qgrams_ = new unsigned char[capacity_ * QGRAM_VECTOR_BYTES];
        memcpy(qgrams_, other.qgrams_, sizeof(unsigned char) * capacity_ * QGRAM_VECTOR_BYTES);
        qgram_offsets_ = other.qgram_offsets_;

        return *this;

    }

    QgramsDacsArray& QgramsDacsArray::operator=(QgramsDacsArray&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] qgrams_;

        // copy / transfer new resources
        capacity_ = other.capacity_;
        size_ = other.size_;

        qgrams_ = other.qgrams_; other.qgrams_ = nullptr;
        qgram_offsets_ = std::move(other.qgram_offsets_);

        return *this;

    }

    QgramsDacsArray* QgramsDacsArray::clone() const { // deep-copy clone method
        return new QgramsDacsArray(*this);
    }

    void QgramsDacsArray::push_back(const std::string& seq) {

        extend();

        size_t offset = size_ * QGRAM_VECTOR_BYTES;
        qgram_offsets_.push_back(offset);
        memset(qgrams_ + offset, 0, QGRAM_VECTOR_BYTES);
        unsigned long qgram = 0;
        unsigned long j = 0;
        while((j < QGRAM_LENGTH - 1) && (j < seq.length())) {

            qgram = (qgram << 2) | (acgtu_map[seq[j]] - 1);
            j++;

        }
        while(j < seq.length()) {

            qgram = (qgram << 2) | (acgtu_map[seq[j]] - 1);
            (qgrams_ + offset)[(qgram >> 3) & (QGRAM_VECTOR_BYTES - 1)] ^= (1 << (qgram & 7));
            j++;

        }

        size_++;

    }

    unsigned char* QgramsDacsArray::get(const numSeqs_t i) const {
        return qgrams_ + qgram_offsets_.get(i);
    }

    void QgramsDacsArray::extend() {

        if (size_ == capacity_) {

            unsigned char* tmp4 = new unsigned char[2 * capacity_ * QGRAM_VECTOR_BYTES];
            memcpy(tmp4, qgrams_, capacity_ * QGRAM_VECTOR_BYTES);
            delete[] qgrams_;
            qgrams_ = tmp4;

            capacity_ *= 2;

        }

    }

    void QgramsDacsArray::finalise() {
        qgram_offsets_.finalise();
    }

    void QgramsDacsArray::sort(const std::vector<numSeqs_t>& permutation) {
        qgram_offsets_.sort(permutation);
    }

    size_t QgramsDacsArray::size_in_bytes() const {
        return sizeof(QgramsDacsArray) // QgramsDacsArray itself
            + sizeof(unsigned char) * capacity_ * QGRAM_VECTOR_BYTES // qgrams_
            + qgram_offsets_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>); // qgram_offsets_ (do not count size of DacsCapacityArray object twice)
    }

    void QgramsDacsArray::show_memory() const {

        std::cout << "#  * QgramsDacsArray " << std::endl;
        std::cout << "#  * " << std::endl;
        std::cout << "#  * sizeof(QgramsDacsArray): " << sizeof(QgramsDacsArray) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char): " << sizeof(unsigned char) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char*): " << sizeof(unsigned char*) << " bytes" << std::endl;
        std::cout << "#  * sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>): " << sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned long long): " << sizeof(unsigned long long) << " bytes" << std::endl;
        std::cout << "#  * " << std::endl;
        std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
        std::cout << "#  * qgrams_: " << (sizeof(unsigned char) * capacity_ * QGRAM_VECTOR_BYTES) << " bytes" << std::endl;
        std::cout << "#  * qgram_offsets_: " << qgram_offsets_.size_in_bytes() << " bytes" << std::endl;

    }

}
