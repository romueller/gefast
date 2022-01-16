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

#ifndef GEFAST_NUMERICARRAYS_HPP
#define GEFAST_NUMERICARRAYS_HPP

#include "../Base.hpp"
#include "../modes/SpaceLevenshteinMode.hpp"
#include "Basics.hpp"
#include "dacs/DacsBuilder.hpp"
#include "dacs/Dacs.hpp"
#include "dacs/SplitDacsBuilder.hpp"
#include "dacs/SplitDacs.hpp"

namespace GeFaST {

    /*
     * Basic implementation of NumericComponent storing the values as a C-style array.
     */
    template<typename V>
    class SimpleArray : public NumericComponent<V> {
    public:
        SimpleArray() = default;

        explicit SimpleArray(const unsigned long long capacity) {

            size_ = 0;
            capacity_ = capacity;
            arr_ = new V[capacity_];

        }

        ~SimpleArray() {
            delete[] arr_;
        }

        SimpleArray(const SimpleArray& other) { // copy constructor

            capacity_ = other.capacity_;
            size_ = other.size_;

            arr_ = new V[capacity_];
            for (auto i = 0; i < size_; i++) {
                arr_[i] = other.arr_[i];
            }

        }

        SimpleArray(SimpleArray&& other) noexcept { // move constructor

            capacity_ = other.capacity_;
            size_ = other.size_;

            arr_ = other.arr_; other.arr_ = nullptr;

        }

        SimpleArray& operator=(const SimpleArray& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] arr_;

            // copy new resources
            capacity_ = other.capacity_;
            size_ = other.size_;

            arr_ = new V[capacity_];
            for (auto i = 0; i < size_; i++) {
                arr_[i] = other.arr_[i];
            }

            return *this;

        }

        SimpleArray& operator=(SimpleArray&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] arr_;

            // copy / transfer new resources
            capacity_ = other.capacity_;
            size_ = other.size_;

            arr_ = other.arr_; other.arr_ = nullptr;

            return *this;

        }

        SimpleArray* clone() const override { // deep-copy clone method
            return new SimpleArray(*this);
        }

        inline void push_back(const V v) override {

            extend();

            arr_[size_] = v;
            size_++;

        }

        inline V get(const numSeqs_t i) const override {
            return arr_[i];
        }

        inline void extend() override {

            if (size_ == capacity_) {

                V* tmp = new V[2 * capacity_];
                memcpy(tmp, arr_, capacity_ * sizeof(V));
                delete[] arr_;
                arr_ = tmp;

                capacity_ *= 2;

            }

        }

        void finalise() override {
            // nothing to do
        }

		void sort(const std::vector<numSeqs_t>& permutation) override {

			std::vector<bool> done(size_);
            for (numSeqs_t i = 0; i < size_; ++i) {

                if (!done[i]) {

                    done[i] = true;
                    numSeqs_t prev_j = i;
                    numSeqs_t j = permutation[i];

                    while (i != j) {

						//swap_values(prev_j, j);
						std::swap(arr_[prev_j], arr_[j]);

                        done[j] = true;
                        prev_j = j;
                        j = permutation[j];

                    }

                }

            }

		}

        size_t size_in_bytes() const override {
            return sizeof(SimpleArray<V>) // SimpleArray itself
                + sizeof(V) * capacity_; // arr_
        }

        void show_memory() const override {

            std::cout << "#  * SimpleArray " << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * sizeof(SimpleArray<V>): " << sizeof(SimpleArray<V>) << " bytes" << std::endl;
            std::cout << "#  * sizeof(V): " << sizeof(V) << " bytes" << std::endl;
            std::cout << "#  * sizeof(unsigned long long): " << sizeof(unsigned long long) << " bytes" << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
            std::cout << "#  * arr_: " << (sizeof(V) * capacity_) << " bytes" << std::endl;

        }

    private:
        V* arr_ = nullptr; // array of stored (numeric) values

        unsigned long long capacity_ = 0; // capacity of the collection (in terms of the number of amplicons)
        unsigned long long size_ = 0; // current number of amplicons in the collection

    };

    /*
     * Basic implementation of QgramComponent for storing the q-gram profiles
     * as bit vectors like in, e.g., ArrayQgramAmpliconCollection.
     */
    class QgramsArray : public QgramComponent {
    public:
        QgramsArray() = default;

        QgramsArray(const unsigned long long capacity);

        ~QgramsArray();

        QgramsArray(const QgramsArray& other);

        QgramsArray(QgramsArray&& other) noexcept;

        QgramsArray& operator=(const QgramsArray& other);

        QgramsArray& operator=(QgramsArray&& other) noexcept;

        QgramsArray* clone() const override;

        void push_back(const std::string& seq) override;

        unsigned char* get(const numSeqs_t i) const override;

        void extend() override;

        void finalise() override;

		void sort(const std::vector<numSeqs_t>& permutation) override;

        size_t size_in_bytes() const override;

        void show_memory() const override;

    private:
        unsigned char* qgrams_ = nullptr; // array of q-gram profiles
        numSeqs_t* qgram_offsets_ = nullptr; // offsets to start positions of the individual q-gram profiles in qgrams_

        unsigned long long capacity_ = 0; // capacity of the collection (in terms of the number of amplicons)
        unsigned long long size_ = 0; // current number of amplicons in the collection

    };



    /*
     * General implementation of NumericComponent storing the values of type V in a DACs sequence of type D.
     * The DACs sequence is constructed using a (temporary) builder of type B.
     * Can be used with Dacs and SplitDacs.
     *
     * Requires to know the distribution of values to be stored to appropriately initialise the builder.
     * Once all values have been added, the DacsArray is finalised by obtaining the actual DACs sequence from the builder.
     */
    template<typename D, typename B, typename V>
    class DacsArray : public NumericComponent<V> {
    public:
        DacsArray() = default;

        explicit DacsArray(const std::map<ularge_t, ularge_t>& counts) {

            // set up builder based on distribution
            builder_ = new B(counts);
            arr_ = nullptr;

        }

        ~DacsArray() {

            delete arr_;
            delete builder_;

        }

        DacsArray(const DacsArray& other) { // copy constructor

            arr_ = (arr_ == nullptr) ? nullptr : other.arr_->clone();
            builder_ = (builder_ == nullptr) ? nullptr : other.builder_->clone();

        }

        DacsArray(DacsArray&& other) noexcept { // move constructor

            arr_ = other.arr_; other.arr_ = nullptr;
            builder_ = other.builder_; other.builder_ = nullptr;

        }

        DacsArray& operator=(const DacsArray& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete arr_;
            delete builder_;

            // copy new resources
            arr_ = (arr_ == nullptr) ? nullptr : other.arr_->clone();
            builder_ = (builder_ == nullptr) ? nullptr : other.builder_->clone();

            return *this;

        }

        DacsArray& operator=(DacsArray&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete arr_;
            delete builder_;

            // copy / transfer new resources
            arr_ = other.arr_; other.arr_ = nullptr;
            builder_ = other.builder_; other.builder_ = nullptr;

            return *this;

        }

        DacsArray* clone() const override { // deep-copy clone method
            return new DacsArray(*this);
        }

        inline void push_back(const V v) override {
            builder_->append(v);
        }

        inline V get(const numSeqs_t i) const override {
            return (*arr_)[i];
        }

        inline void extend() override {
            // nothing to do
        }

        void finalise() override {

            if (builder_ != nullptr) {

                arr_ = builder_->finalise();
                delete builder_; builder_ = nullptr;

            }

        }

        void sort(const std::vector<numSeqs_t>& permutation) override {

            // build actual DACs instance of unsorted values
            finalise();

            // build sorted DACs instance and replace unsorted one
            auto tmp = new D(*arr_, permutation);
            delete arr_;
            arr_ = tmp;

        }

        size_t size_in_bytes() const override {
            return sizeof(DacsArray<D, B, V>) // DacsArray itself
                + ((builder_ != nullptr) ? builder_->size_in_bytes() : 0) // builder_
                + ((arr_ != nullptr) ? arr_->size_in_bytes() : 0); // arr_
        }

        void show_memory() const override {

            std::cout << "#  * DacsArray " << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * sizeof(DacsArray<D, B, V>): " << sizeof(DacsArray<D, B, V>) << " bytes" << std::endl;
            std::cout << "#  * sizeof(D*): " << sizeof(D*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(B*): " << sizeof(B*) << " bytes" << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
            std::cout << "#  * builder_: " << ((builder_ != nullptr) ? builder_->size_in_bytes() : 0) << " bytes" << std::endl;
            std::cout << "#  * arr_: " << ((arr_ != nullptr) ? arr_->size_in_bytes() : 0) << " bytes" << std::endl;

        }

    private:
        B* builder_ = nullptr; // (temporary) builder for DACs sequence
        D* arr_ = nullptr; // array of stored numeric values

    };

    /*
     * Specialised implementation of NumericComponent storing the values of type V in a DACs sequence of type D.
     * The DACs sequence is constructed using a (temporary) builder of type B.
     * Can be used with Dacs and SplitDacs.
     *
     * Designed for cases where the minimum value is quite large.
     * Internally, the minimum is subtracted from each value before inserting it into the DACs sequence
     * and added again during retrieval.
     *
     * Requires to know the distribution of values to be stored to appropriately initialise the builder.
     * Once all values have been added, the DacsArray is finalised by obtaining the actual DACs sequence from the builder.
     */
    template<typename D, typename B, typename V>
    class DacsOffsetArray : public NumericComponent<V> {
    public:
        DacsOffsetArray() = default;

        explicit DacsOffsetArray(const std::map<ularge_t, ularge_t>& counts) {

            // reduce all length values by offset
            std::map<ularge_t, ularge_t> changed;
            offset_ = counts.begin()->first; // = shortest length
            for (auto& kv : counts) {
                changed[kv.first - offset_] = kv.second;
            }

            // set up builder based on distribution
            builder_ = new B(changed);
            arr_ = nullptr;

        }

        ~DacsOffsetArray() {

            delete arr_;
            delete builder_;

        }

        DacsOffsetArray(const DacsOffsetArray& other) { // copy constructor

            offset_ = other.offset_;

            arr_ = (arr_ == nullptr) ? nullptr : other.arr_->clone();
            builder_ = (builder_ == nullptr) ? nullptr : other.builder_->clone();

        }

        DacsOffsetArray(DacsOffsetArray&& other) noexcept { // move constructor

            offset_ = other.offset_;

            arr_ = other.arr_; other.arr_ = nullptr;
            builder_ = other.builder_; other.builder_ = nullptr;

        }

        DacsOffsetArray& operator=(const DacsOffsetArray& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete arr_;
            delete builder_;

            // copy new resources
            offset_ = other.offset_;

            arr_ = (arr_ == nullptr) ? nullptr : other.arr_->clone();
            builder_ = (builder_ == nullptr) ? nullptr : other.builder_->clone();

            return *this;

        }

        DacsOffsetArray& operator=(DacsOffsetArray&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete arr_;
            delete builder_;

            // copy / transfer new resources
            offset_ = other.offset_;

            arr_ = other.arr_; other.arr_ = nullptr;
            builder_ = other.builder_; other.builder_ = nullptr;

            return *this;

        }

        DacsOffsetArray* clone() const override { // deep-copy clone method
            return new DacsOffsetArray(*this);
        }

        inline void push_back(const V v) override {
            builder_->append(v - offset_);
        }

        inline V get(const numSeqs_t i) const override {
            return (*arr_)[i] + offset_;
        }

        inline void extend() override {
            // nothing to do
        }

        void finalise() override {

            if (builder_ != nullptr) {

                arr_ = builder_->finalise();
                delete builder_; builder_ = nullptr;

            }

        }

        void sort(const std::vector<numSeqs_t>& permutation) override {

            // build actual DACs instance of unsorted values
            finalise();

            // build sorted DACs instance and replace unsorted one
            auto tmp = new D(*arr_, permutation);
            delete arr_;
            arr_ = tmp;

        }

        size_t size_in_bytes() const override {
            return sizeof(DacsOffsetArray<D, B, V>) // DacsOffsetArray itself
                + ((builder_ != nullptr) ? builder_->size_in_bytes() : 0) // builder_
                + ((arr_ != nullptr) ? arr_->size_in_bytes() : 0); // arr_
        }

        void show_memory() const override {

            std::cout << "#  * DacsOffsetArray " << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * sizeof(DacsOffsetArray<D, B, V>): " << sizeof(DacsArray<D, B, V>) << " bytes" << std::endl;
            std::cout << "#  * sizeof(D*): " << sizeof(D*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(B*): " << sizeof(B*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(lenSeqs_t): " << sizeof(lenSeqs_t) << " bytes" << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
            std::cout << "#  * builder_: " << ((builder_ != nullptr) ? builder_->size_in_bytes() : 0) << " bytes" << std::endl;
            std::cout << "#  * arr_: " << ((arr_ != nullptr) ? arr_->size_in_bytes() : 0) << " bytes" << std::endl;

        }

    private:
        B* builder_ = nullptr; // (temporary) builder for DACs sequence
        D* arr_ = nullptr; // array of stored numeric values
        lenSeqs_t offset_ = 0; // offset subtracted from (added to) all length values before insertion (after retrieval)

    };

    /*
     * Specialised implementation of NumericComponent storing the values of type V in a DACs sequence of type D.
     * The DACs sequence is constructed using a (temporary) builder of type B.
     * Can be used with Dacs and SplitDacs.
     *
     * Designed for cases where a large proportion of the values is 1.
     * Internally, only values larger than 1 are stored in the DACs sequence.
     * In addition, a bit vector indicates whether a retrieval should look into the DACs sequence or directly return 1.
     * Uses a rank data structure to map positions in the overall sequence to the DACs sequence.
     *
     * Requires to know the distribution of values to be stored to appropriately initialise the builder.
     * Once all values have been added, the DacsArray is finalised by obtaining the actual DACs sequence from the builder.
     */
    template<typename D, typename B, typename V, typename S>
    class DacsBitsArray : public NumericComponent<V> {
    public:
        DacsBitsArray() = default;

        explicit DacsBitsArray(const std::map<ularge_t, ularge_t>& counts) : ranks_(BitRank<S, 20>(nullptr, 0, false)) {

            // determine distribution of values > 1 for initialisation of builder
            std::map<ularge_t, ularge_t> changed;
            size_ = 0;
            for (auto& kv : counts) {

                if (kv.first != 1) {
                    changed[kv.first] = kv.second;
                }
                size_ += kv.second;

            }

            // set up flag data structure and set size_ to current size
            flags_ = BitVector<S>(size_);
            size_ = 0;

            // set up builder and DACs depending on the existence of at least one value > 1
            builder_ = (changed.size() > 0) ? new B(changed) : nullptr;
            arr_ = (changed.size() > 0) ? nullptr : new D();

        }

        ~DacsBitsArray() {

            delete arr_;
            delete builder_;

        }

        DacsBitsArray(const DacsBitsArray& other) : flags_(other.flags_), ranks_(flags_) { // copy constructor

            size_ = other.size_;

            arr_ = (arr_ == nullptr) ? nullptr : other.arr_->clone();
            builder_ = (builder_ == nullptr) ? nullptr : other.builder_->clone();

        }

        DacsBitsArray(DacsBitsArray&& other) noexcept : flags_(std::move(other.flags_)), ranks_(flags_) { // move constructor

            size_ = other.size_;

            arr_ = other.arr_; other.arr_ = nullptr;
            builder_ = other.builder_; other.builder_ = nullptr;

        }

        DacsBitsArray& operator=(const DacsBitsArray& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete arr_;
            delete builder_;

            // copy new resources
            size_ = other.size_;

            arr_ = (arr_ == nullptr) ? nullptr : other.arr_->clone();
            builder_ = (builder_ == nullptr) ? nullptr : other.builder_->clone();
            flags_ = other.flags_;
            ranks_ = BitRank<S, 20>(flags_);

            return *this;

        }

        DacsBitsArray& operator=(DacsBitsArray&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete arr_;
            delete builder_;

            // copy / transfer new resources
            size_ = other.size_;

            arr_ = other.arr_; other.arr_ = nullptr;
            builder_ = other.builder_; other.builder_ = nullptr;
            flags_ = other.flags_;
            ranks_ = BitRank<S, 20>(flags_);

            return *this;

        }

        DacsBitsArray* clone() const override { // deep-copy clone method
            return new DacsBitsArray(*this);
        }

        inline void push_back(const V v) override {

            if (v > 1) {
                builder_->append(v);
            }

            flags_.set(size_, v > 1);
            size_++;

        }

        inline V get(const numSeqs_t i) const override {
            return (flags_[i]) ? (*arr_)[rank(i) - 1] : 1;
        }

        inline void extend() override {
            // nothing to do
        }

        void finalise() override {

            if (builder_ != nullptr) {

                arr_ = builder_->finalise();
                delete builder_; builder_ = nullptr;

                ranks_ = BitRank<S, 20>(flags_);

            }

        }

        void sort(const std::vector<numSeqs_t>& permutation) override {

            // build actual DACs instance (effect only if not empty)
            finalise();

            // temporarily determine distribution of values > 1 for builder of sorted sequence
            {
                std::map<ularge_t, ularge_t> changed;
                auto reader = arr_->reader();
                while (reader.has_next()) {
                    changed[reader.next()]++;
                }
                builder_ = (changed.size() > 0) ? new B(changed) : nullptr;
            }

            // if at least one value > 1, build the sorted DACs (otherwise nothing has to be done)
            if (builder_ != nullptr) {

                std::vector<bool> tmp_flags;

                for (auto i = 0; i < permutation.size(); i++) {

                    V val = get(permutation[i]);
                    if (val > 1) {
                        builder_->append(val);
                    }
                    tmp_flags.push_back(val > 1);

                }

                delete arr_;
                flags_.set(0, tmp_flags);

            }

            // build actual sorted DACs instance (effect only if not empty)
            finalise();

        }

        size_t size_in_bytes() const override {
            return sizeof(DacsBitsArray<D, B, V, S>) // DacsBitsArray itself
                + ((builder_ != nullptr) ? builder_->size_in_bytes() : 0) // builder_
                + ((arr_ != nullptr) ? arr_->size_in_bytes() : 0) // arr_
                + flags_.size_in_bytes() - sizeof(BitVector<S>) // flags_ (do not count size of BitVector object twice)
                + ranks_.size_in_bytes() - sizeof(BitRank<S, 20>);  // ranks_ (do not count size of BitRank object twice)
        }

        void show_memory() const override {

            std::cout << "#  * DacsBitsArray " << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * sizeof(DacsBitsArray<D, B, V, S>): " << sizeof(DacsArray<D, B, V>) << " bytes" << std::endl;
            std::cout << "#  * sizeof(D*): " << sizeof(D*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(B*): " << sizeof(B*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(numSeqs_t): " << sizeof(numSeqs_t) << " bytes" << std::endl;
            std::cout << "#  * sizeof(BitVector<S>): " << sizeof(BitVector<S>) << " bytes" << std::endl;
            std::cout << "#  * sizeof(BitRank<S, 20>): " << sizeof(BitRank<S, 20>) << " bytes" << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
            std::cout << "#  * builder_: " << ((builder_ != nullptr) ? builder_->size_in_bytes() : 0) << " bytes" << std::endl;
            std::cout << "#  * arr_: " << ((arr_ != nullptr) ? arr_->size_in_bytes() : 0) << " bytes" << std::endl;
            std::cout << "#  * flags_: " << flags_.size_in_bytes() << " bytes" << std::endl;
            std::cout << "#  * ranks_: " << ranks_.size_in_bytes() << " bytes" << std::endl;

        }

    private:
        B* builder_ = nullptr; // (temporary) builder for DACs sequence
        D* arr_ = nullptr; // array of stored numeric values
        numSeqs_t size_ = 0; // current size of the array
        BitVector<S> flags_; // flags indicating whether the value is > 1 (stored in DACs sequence)
        BitRank<S, 20> ranks_; // rank data structure on top of flags_

        // determine rank among values > 1 (to access DACs sequence)
        numSeqs_t rank(numSeqs_t i) const {
            return ranks_.rank(i);
        }

    };



    /*
     * General implementation of NumericComponent storing the values of type V in a DACs sequence of type D.
     * The DACs sequence is constructed using a (temporary) builder of type B.
     * Can only be used with SplitDacs.
     *
     * Does not require to know the distribution of values to be stored.
     * The builder uses an initial capacity and extends the DACs representation as needed
     * by enlarging levels ((new size) = (current size) * (extend_factor))
     * and adding further levels ((size of level) = (size of previous level) * (level_factor)).
     * Uses the same chunk length on each level but during the construction of the actual DACs instance
     * an optional optimisation can be performed.
     * Once all values have been added, the DacsArray is finalised by obtaining the actual DACs sequence from the builder.
     */
    template<typename D, typename B, typename V>
    class DacsCapacityArray : public NumericComponent<V> {
    public:
        DacsCapacityArray() = default;

        explicit DacsCapacityArray(usmall_t chunk_length, ularge_t capacity, float level_factor, float extend_factor, bool opt) {

            // set up builder based on distribution
            builder_ = new B(chunk_length, capacity, level_factor, extend_factor);
            arr_ = nullptr;
            opt_ = opt;

        }

        ~DacsCapacityArray() {

            delete arr_;
            delete builder_;

        }

        DacsCapacityArray(const DacsCapacityArray& other) { // copy constructor

            arr_ = (arr_ == nullptr) ? nullptr : other.arr_->clone();
            builder_ = (builder_ == nullptr) ? nullptr : other.builder_->clone();
            opt_ = other.opt_;

        }

        DacsCapacityArray(DacsCapacityArray&& other) noexcept { // move constructor

            arr_ = other.arr_; other.arr_ = nullptr;
            builder_ = other.builder_; other.builder_ = nullptr;
            opt_ = other.opt_;

        }

        DacsCapacityArray& operator=(const DacsCapacityArray& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete arr_;
            delete builder_;

            // copy new resources
            arr_ = (arr_ == nullptr) ? nullptr : other.arr_->clone();
            builder_ = (builder_ == nullptr) ? nullptr : other.builder_->clone();
            opt_ = other.opt_;

            return *this;

        }

        DacsCapacityArray& operator=(DacsCapacityArray&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete arr_;
            delete builder_;

            // copy / transfer new resources
            arr_ = other.arr_; other.arr_ = nullptr;
            builder_ = other.builder_; other.builder_ = nullptr;
            opt_ = other.opt_;

            return *this;

        }

        DacsCapacityArray* clone() const override { // deep-copy clone method
            return new DacsCapacityArray(*this);
        }

        inline void push_back(const V v) override {
            builder_->append(v);
        }

        inline V get(const numSeqs_t i) const override {
            return (*arr_)[i];
        }

        inline void extend() override {
            // nothing to do
        }

        void finalise() override {

            if (builder_ != nullptr) {

                arr_ = builder_->finalise(opt_);
                delete builder_; builder_ = nullptr;

            }

        }

        bool finalised() const {
            return builder_ == nullptr;
        }

        void sort(const std::vector<numSeqs_t>& permutation) override {

            // build actual DACs instance of unsorted values
            finalise();

            // build sorted DACs instance and replace unsorted one
            auto tmp = new D(*arr_, permutation);
            delete arr_;
            arr_ = tmp;

        }

        size_t size_in_bytes() const override {
            return sizeof(DacsCapacityArray<D, B, V>) // DacsCapacityArray itself
                + ((builder_ != nullptr) ? builder_->size_in_bytes() : 0) // builder_
                + ((arr_ != nullptr) ? arr_->size_in_bytes() : 0); // arr_
        }

        void show_memory() const override {

            std::cout << "#  * DacsCapacityArray " << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * sizeof(DacsCapacityArray<D, B, V>): " << sizeof(DacsArray<D, B, V>) << " bytes" << std::endl;
            std::cout << "#  * sizeof(D*): " << sizeof(D*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(B*): " << sizeof(B*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(bool): " << sizeof(bool) << " bytes" << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
            std::cout << "#  * builder_: " << ((builder_ != nullptr) ? builder_->size_in_bytes() : 0) << " bytes" << std::endl;
            std::cout << "#  * arr_: " << ((arr_ != nullptr) ? arr_->size_in_bytes() : 0) << " bytes" << std::endl;

        }

    private:
        B* builder_ = nullptr; // (temporary) builder for DACs sequence
        D* arr_ = nullptr; // array of stored numeric values
        bool opt_ = false; // flag indicating whether the DACs instance should be optimised during the finalisation

    };




    /*
     * Variant of QgramsArray: the offsets to the different profiles are stored as a DACs sequence.
     */
    class QgramsDacsArray : public QgramComponent {
    public:
        QgramsDacsArray() = default;

        QgramsDacsArray(const unsigned long long capacity, const ChunkLengths& config);

        ~QgramsDacsArray();

        QgramsDacsArray(const QgramsDacsArray& other);

        QgramsDacsArray(QgramsDacsArray&& other) noexcept;

        QgramsDacsArray& operator=(const QgramsDacsArray& other);

        QgramsDacsArray& operator=(QgramsDacsArray&& other) noexcept;

        QgramsDacsArray* clone() const override;

        void push_back(const std::string& seq) override;

        unsigned char* get(const numSeqs_t i) const override;

        void extend() override;

        void finalise() override;

        void sort(const std::vector<numSeqs_t>& permutation) override;

        size_t size_in_bytes() const override;

        void show_memory() const override;

    private:
        typedef SplitDacs<numSeqs_t, size_t> dacs_t;
        typedef SplitDacsDirectBuilder<numSeqs_t, size_t> builder_t;

        unsigned char* qgrams_ = nullptr; // array of q-gram profiles
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> qgram_offsets_; // offsets to start positions of the individual q-gram profiles in qgrams_

        unsigned long long capacity_ = 0; // capacity of the collection (in terms of the number of amplicons)
        unsigned long long size_ = 0; // current number of amplicons in the collection

    };

}

#endif //GEFAST_NUMERICARRAYS_HPP
