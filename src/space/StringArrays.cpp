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

#include "../../include/space/StringArrays.hpp"

namespace GeFaST {

    /* === CollectiveStringArray === */

    CollectiveStringArray::CollectiveStringArray(const unsigned long long capacity, const unsigned long long capacity_strings) {

        size_ = 0;
        capacity_ = capacity + 1; // + 1 to have space for a sentinel after the currently last entry
        capacity_strings_ = capacity_strings + 1; // + 1 to have space for a sentinel after the currently last entry

        strings_ = new char[capacity_strings_];
        string_pointers_ = new char*[capacity_];
        next_string_ = strings_;

    }

    CollectiveStringArray::~CollectiveStringArray() {

        delete[] string_pointers_;
        delete[] strings_;

    }

    CollectiveStringArray::CollectiveStringArray(const CollectiveStringArray& other) { // copy constructor

        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_strings_ = other.capacity_strings_;
        strings_ = new char[capacity_strings_];
        memcpy(strings_, other.strings_, sizeof(char) * capacity_strings_);
        string_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            string_pointers_[i] = strings_ + (other.string_pointers_[i] - other.strings_);
        }
        next_string_ = strings_ + (other.next_string_ - other.strings_);

    }

    CollectiveStringArray::CollectiveStringArray(CollectiveStringArray&& other) noexcept { // move constructor

        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_strings_ = other.capacity_strings_;
        next_string_ = other.next_string_;
        strings_ = other.strings_; other.strings_ = nullptr;
        string_pointers_ = other.string_pointers_; other.string_pointers_ = nullptr;

    }

    CollectiveStringArray& CollectiveStringArray::operator=(const CollectiveStringArray& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] string_pointers_;
        delete[] strings_;

        // copy new resources
        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_strings_ = other.capacity_strings_;
        strings_ = new char[capacity_strings_];
        memcpy(strings_, other.strings_, sizeof(char) * capacity_strings_);
        string_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            string_pointers_[i] = strings_ + (other.string_pointers_[i] - other.strings_);
        }
        next_string_ = strings_ + (other.next_string_ - other.strings_);

        return *this;

    }

    CollectiveStringArray& CollectiveStringArray::operator=(CollectiveStringArray&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] string_pointers_;
        delete[] strings_;

        // copy / transfer new resources
        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_strings_ = other.capacity_strings_;
        next_string_ = other.next_string_;
        strings_ = other.strings_; other.strings_ = nullptr;
        string_pointers_ = other.string_pointers_; other.string_pointers_ = nullptr;

        return *this;

    }

    CollectiveStringArray* CollectiveStringArray::clone() const { // deep-copy clone method
        return new CollectiveStringArray(*this);
    }

    void CollectiveStringArray::push_back(const std::string& s) {

        extend(s);

        string_pointers_[size_] = next_string_;
        strcpy(next_string_, s.c_str());
        next_string_ += s.length() + 1;

        size_++;
        string_pointers_[size_] = next_string_; // sentinel for length calculations

    }

    std::unique_ptr<SequenceWrapper> CollectiveStringArray::get(const numSeqs_t i, const lenSeqs_t len) const {
        return std::unique_ptr<SequenceWrapper>(new PointerSequenceWrapper(string_pointers_[i], len));
    }

    std::string CollectiveStringArray::get_str(const numSeqs_t i, const lenSeqs_t len) const {
        return std::string(string_pointers_[i]);
    }

    void CollectiveStringArray::extend(const std::string& s) {

        size_t l = s.length();
        unsigned long long len = next_string_ - strings_;
        if (len + l + 1 > capacity_strings_) {

            char* tmp1 = new char[2 * capacity_strings_];
            memcpy(tmp1, strings_, len * sizeof(char));
            char** tmp2 = new char*[capacity_];
            for (auto i = 0; i < size_; i++) {
                tmp2[i] = tmp1 + (string_pointers_[i] - strings_);
            }
            delete[] strings_;
            delete[] string_pointers_;
            strings_ = tmp1;
            next_string_ = strings_ + len;
            string_pointers_ = tmp2;

            capacity_strings_ = std::max(len + l + 1, 2 * capacity_strings_);

        }

        if (size_ + 1 == capacity_) {

            char** tmp1 = new char*[2 * capacity_];
            memcpy(tmp1, string_pointers_, capacity_ * sizeof(char*));
            delete[] string_pointers_;
            string_pointers_ = tmp1;

            capacity_ *= 2;

        }

    }

    void CollectiveStringArray::finalise() {
        // nothing to do
    }

    void CollectiveStringArray::sort(const std::vector<numSeqs_t>& permutation) {

        std::vector<bool> done(size_);
        for (numSeqs_t i = 0; i < size_; ++i) {

            if (!done[i]) {

                done[i] = true;
                numSeqs_t prev_j = i;
                numSeqs_t j = permutation[i];

                while (i != j) {

                    // swap_values(prev_j, j);
                    std::swap(string_pointers_[prev_j], string_pointers_[j]);

                    done[j] = true;
                    prev_j = j;
                    j = permutation[j];

                }

            }

        }

    }

    size_t CollectiveStringArray::size_in_bytes() const {
        return sizeof(CollectiveStringArray) // CollectiveStringArray itself
            + sizeof(char) * capacity_strings_ // strings_
            + sizeof(char*) * capacity_; // string_pointers_
    }

    void CollectiveStringArray::show_memory() const {

        std::cout << "#  * CollectiveStringArray " << std::endl;
        std::cout << "#  * " << std::endl;
        std::cout << "#  * sizeof(CollectiveStringArray): " << sizeof(CollectiveStringArray) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char): " << sizeof(unsigned char) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char*): " << sizeof(unsigned char*) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char**): " << sizeof(unsigned char**) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned long long): " << sizeof(unsigned long long) << " bytes" << std::endl;
        std::cout << "#  * sizeof(numSeqs_t): " << sizeof(numSeqs_t) << " bytes" << std::endl;
        std::cout << "#  * " << std::endl;
        std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
        std::cout << "#  * strings_: " << (sizeof(char) * capacity_strings_) << " bytes" << std::endl;
        std::cout << "#  * string_pointers_: " << (sizeof(char*) * capacity_) << " bytes" << std::endl;

    }



    /* === PrefixStringArray === */

    PrefixStringArray::PrefixStringArray(const unsigned char p, const unsigned long long capacity, const unsigned long long capacity_strings) {

        size_ = 0;
        capacity_ = capacity + 1; // + 1 to have space for a sentinel after the currently last entry
        capacity_strings_ = capacity_strings - capacity * p + 1; // + 1 to have space for a sentinel after the currently last entry

        strings_ = new char[capacity_strings_];
        string_pointers_ = new char*[capacity_];
        next_string_ = strings_;

        prefix_len_ = p;
        prefixes_ = new unsigned char[capacity_];

    }

    PrefixStringArray::~PrefixStringArray() {

        delete[] string_pointers_;
        delete[] strings_;
        delete[] prefixes_;

    }

    PrefixStringArray::PrefixStringArray(const PrefixStringArray& other) { // copy constructor

        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_strings_ = other.capacity_strings_;
        strings_ = new char[capacity_strings_];
        memcpy(strings_, other.strings_, sizeof(char) * capacity_strings_);
        string_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            string_pointers_[i] = strings_ + (other.string_pointers_[i] - other.strings_);
        }
        next_string_ = strings_ + (other.next_string_ - other.strings_);

        prefix_len_ = other.prefix_len_;
        pref_num_map_ = other.pref_num_map_;
        num_pref_map_ = other.num_pref_map_;
        prefixes_ = new unsigned char[capacity_];
        for (auto i = 0; i < size_; i++) {
            prefixes_[i] = other.prefixes_[i];
        }

    }

    PrefixStringArray::PrefixStringArray(PrefixStringArray&& other) noexcept { // move constructor

        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_strings_ = other.capacity_strings_;
        next_string_ = other.next_string_;
        strings_ = other.strings_; other.strings_ = nullptr;
        string_pointers_ = other.string_pointers_; other.string_pointers_ = nullptr;

        prefix_len_ = other.prefix_len_;
        pref_num_map_ = other.pref_num_map_;
        num_pref_map_ = other.num_pref_map_;
        prefixes_ = other.prefixes_; other.prefixes_ = nullptr;

    }

    PrefixStringArray& PrefixStringArray::operator=(const PrefixStringArray& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] string_pointers_;
        delete[] strings_;
        delete[] prefixes_;

        // copy new resources
        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_strings_ = other.capacity_strings_;
        strings_ = new char[capacity_strings_];
        memcpy(strings_, other.strings_, sizeof(char) * capacity_strings_);
        string_pointers_ = new char*[capacity_];
        for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
            string_pointers_[i] = strings_ + (other.string_pointers_[i] - other.strings_);
        }
        next_string_ = strings_ + (other.next_string_ - other.strings_);

        prefix_len_ = other.prefix_len_;
        pref_num_map_ = other.pref_num_map_;
        num_pref_map_ = other.num_pref_map_;
        prefixes_ = new unsigned char[capacity_];
        for (auto i = 0; i < size_; i++) {
            prefixes_[i] = other.prefixes_[i];
        }

        return *this;

    }

    PrefixStringArray& PrefixStringArray::operator=(PrefixStringArray&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] string_pointers_;
        delete[] strings_;
        delete[] prefixes_;

        // copy / transfer new resources
        capacity_ = other.capacity_;
        size_ = other.size_;

        capacity_strings_ = other.capacity_strings_;
        next_string_ = other.next_string_;
        strings_ = other.strings_; other.strings_ = nullptr;
        string_pointers_ = other.string_pointers_; other.string_pointers_ = nullptr;

        prefix_len_ = other.prefix_len_;
        pref_num_map_ = other.pref_num_map_;
        num_pref_map_ = other.num_pref_map_;
        prefixes_ = other.prefixes_; other.prefixes_ = nullptr;

        return *this;

    }

    PrefixStringArray* PrefixStringArray::clone() const { // deep-copy clone method
        return new PrefixStringArray(*this);
    }

    void PrefixStringArray::push_back(const std::string& s) {

        extend(s);

        string_pointers_[size_] = next_string_;
        strcpy(next_string_, s.c_str() + prefix_len_);
        next_string_ += s.length() - prefix_len_ + 1;

        std::string prefix = s.substr(0, prefix_len_);
        auto iter = pref_num_map_.find(prefix);

        if (iter == pref_num_map_.end()) { // previously unseen prefix -> extend mapping

            pref_num_map_[prefix] = num_pref_map_.size();
            prefixes_[size_] = num_pref_map_.size();
            num_pref_map_.push_back(prefix);

        } else { // previously seen prefix -> reuse mapping
            prefixes_[size_] = iter->second;
        }

        size_++;
        string_pointers_[size_] = next_string_; // sentinel for length calculations

    }

    std::unique_ptr<SequenceWrapper> PrefixStringArray::get(const numSeqs_t i, const lenSeqs_t len) const {
        return std::unique_ptr<SequenceWrapper>(new PrefixSeqWrapper(this, i, num_pref_map_[prefixes_[i]].length() + strlen(string_pointers_[i])));
    }

    std::string PrefixStringArray::get_str(const numSeqs_t i, const lenSeqs_t len) const {
        return num_pref_map_[prefixes_[i]] + std::string(string_pointers_[i]);
    }

    void PrefixStringArray::extend(const std::string& s) {

        size_t l = s.length();
        unsigned long long len = next_string_ - strings_;
        if (len + l - prefix_len_ + 1 > capacity_strings_) {

            char* tmp1 = new char[2 * capacity_strings_];
            memcpy(tmp1, strings_, len * sizeof(char));
            char** tmp2 = new char*[capacity_];
            for (auto i = 0; i < size_; i++) {
                tmp2[i] = tmp1 + (string_pointers_[i] - strings_);
            }
            delete[] strings_;
            delete[] string_pointers_;
            strings_ = tmp1;
            next_string_ = strings_ + len;
            string_pointers_ = tmp2;

            capacity_strings_ = std::max(len + l - prefix_len_ + 1, 2 * capacity_strings_);

        }

        if (size_ + 1 == capacity_) {

            char** tmp1 = new char*[2 * capacity_];
            memcpy(tmp1, string_pointers_, capacity_ * sizeof(char*));
            delete[] string_pointers_;
            string_pointers_ = tmp1;

            unsigned char* tmp2 = new unsigned char[2 * capacity_];
            memcpy(tmp2, prefixes_, capacity_ * sizeof(unsigned char));
            delete[] prefixes_;
            prefixes_ = tmp2;

            capacity_ *= 2;

        }

    }

    void PrefixStringArray::finalise() {
        // nothing to do
    }

    void PrefixStringArray::sort(const std::vector<numSeqs_t>& permutation) {

        std::vector<bool> done(size_);
        for (numSeqs_t i = 0; i < size_; ++i) {

            if (!done[i]) {

                done[i] = true;
                numSeqs_t prev_j = i;
                numSeqs_t j = permutation[i];

                while (i != j) {

                    // swap_values(prev_j, j);
                    std::swap(string_pointers_[prev_j], string_pointers_[j]);
                    std::swap(prefixes_[prev_j], prefixes_[j]);

                    done[j] = true;
                    prev_j = j;
                    j = permutation[j];

                }

            }

        }

    }

    size_t PrefixStringArray::size_in_bytes() const {
        return sizeof(PrefixStringArray) // PrefixStringArray itself
            + sizeof(char) * capacity_strings_ // strings_
            + sizeof(char*) * capacity_ // string_pointers_
            + sizeof(unsigned char) * capacity_ // prefixes_
            + map_str_uchar_size_in_bytes(pref_num_map_) - sizeof(std::map<std::string, unsigned char>) // pref_num_map_ (do not count size of map object twice)
            + vec_str_size_in_bytes(num_pref_map_) - sizeof(std::vector<std::string>); // num_pref_map_ (do not count size of vector object twice)
    }

    void PrefixStringArray::show_memory() const {

        std::cout << "#  * PrefixStringArray " << std::endl;
        std::cout << "#  * " << std::endl;
        std::cout << "#  * sizeof(PrefixStringArray): " << sizeof(PrefixStringArray) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char): " << sizeof(unsigned char) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char*): " << sizeof(unsigned char*) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char**): " << sizeof(unsigned char**) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned long long): " << sizeof(unsigned long long) << " bytes" << std::endl;
        std::cout << "#  * sizeof(numSeqs_t): " << sizeof(numSeqs_t) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char): " << sizeof(unsigned char) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char*): " << sizeof(unsigned char*) << " bytes" << std::endl;
        std::cout << "#  * sizeof(std::map<std::string, unsigned char>): " << sizeof(std::map<std::string, unsigned char>) << " bytes" << std::endl;
        std::cout << "#  * sizeof(std::vector<std::string>): " << sizeof(std::vector<std::string>) << " bytes" << std::endl;
        std::cout << "#  * " << std::endl;
        std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
        std::cout << "#  * strings_: " << (sizeof(char) * capacity_strings_) << " bytes" << std::endl;
        std::cout << "#  * string_pointers_: " << (sizeof(char*) * capacity_) << " bytes" << std::endl;
        std::cout << "#  * prefixes_: " << (sizeof(unsigned char) * capacity_) << " bytes" << std::endl;
        std::cout << "#  * pref_num_map_: " << map_str_uchar_size_in_bytes(pref_num_map_) << " bytes" << std::endl;
        std::cout << "#  * num_pref_map_: " << vec_str_size_in_bytes(num_pref_map_) << " bytes" << std::endl;

    }

    /* === PrefixStringArray::PrefixSeqWrapper === */

    PrefixStringArray::PrefixSeqWrapper::PrefixSeqWrapper(const PrefixStringArray* arr, size_t i, lenSeqs_t len) : arr_(arr), index_(i), len_(len) {
        // nothing else to do
    }

    char PrefixStringArray::PrefixSeqWrapper::operator[](size_t i) const {
        return arr_->get_char(index_, i);
    }

    lenSeqs_t PrefixStringArray::PrefixSeqWrapper::length() const {
        return len_;
    }

    std::string PrefixStringArray::PrefixSeqWrapper::to_string() const {
        return arr_->get_prefix(index_) + arr_->get_suffix(index_);
    }

    void PrefixStringArray::PrefixSeqWrapper::write(std::ostream& os) const {
        os << arr_->get_prefix(index_) << arr_->get_suffix(index_);
    }

    PrefixStringArray::PrefixSeqWrapper* PrefixStringArray::PrefixSeqWrapper::clone() const {
        return new PrefixSeqWrapper(*this);
    }

    std::unique_ptr<SubstringWrapper> PrefixStringArray::PrefixSeqWrapper::substr(lenSeqs_t start, lenSeqs_t len) const {
        return std::unique_ptr<SubstringWrapper>(new PrefixSubWrapper(arr_, index_, start, len));
    }

    size_t PrefixStringArray::PrefixSeqWrapper::size_in_bytes() const {
        return sizeof(PrefixSeqWrapper);
    }

    bool PrefixStringArray::PrefixSeqWrapper::less(const SequenceWrapper& other) const {

        auto& o = static_cast<const PrefixSeqWrapper&>(other);

        return (arr_->get_prefix(index_) < o.arr_->get_prefix(o.index_))
                || ((arr_->prefixes_[index_] == o.arr_->prefixes_[o.index_])
                    && (strcmp(arr_->string_pointers_[index_], o.arr_->string_pointers_[o.index_]) < 0));

    }

    bool PrefixStringArray::PrefixSeqWrapper::equals(const SequenceWrapper& other) const {

        auto& o = static_cast<const PrefixSeqWrapper&>(other);

        return (arr_->prefixes_[index_] == o.arr_->prefixes_[o.index_]) &&
            (strcmp(arr_->string_pointers_[index_], o.arr_->string_pointers_[o.index_]) == 0);

    }

    /* === PrefixStringArray::PrefixSubWrapper === */

    PrefixStringArray::PrefixSubWrapper::PrefixSubWrapper(const PrefixStringArray* arr, size_t i, size_t start, lenSeqs_t len)
        : arr_(arr), index_(i), start_(start), len_(len) {
        // nothing else to do
    }

    char PrefixStringArray::PrefixSubWrapper::operator[](size_t i) const {
        return arr_->get_char(index_, start_ + i);
    }

    lenSeqs_t PrefixStringArray::PrefixSubWrapper::length() const {
        return len_;
    }

    std::string PrefixStringArray::PrefixSubWrapper::to_string() const {
        return arr_->get_string(index_, start_, len_);
    }

    PrefixStringArray::PrefixSubWrapper* PrefixStringArray::PrefixSubWrapper::clone() const {
        return new PrefixSubWrapper(*this);
    }

    void PrefixStringArray::PrefixSubWrapper::shift() {
        start_++;
    }

    size_t PrefixStringArray::PrefixSubWrapper::size_in_bytes() const {
        return sizeof(PrefixSubWrapper);
    }

    bool PrefixStringArray::PrefixSubWrapper::equals(const SubstringWrapper& other) const {

        auto& o = static_cast<const PrefixSubWrapper&>(other);

        bool eq = true;
        for (auto i = 0; eq && i < len_; i++) {
            eq = arr_->get_char(index_, start_ + i) == o[i];
        }

        return eq;

    }



    /* === PrefixDacsStringArray === */

    PrefixDacsStringArray::PrefixDacsStringArray(const unsigned char p, const unsigned long long capacity,
                                                 const unsigned long long capacity_strings, const ChunkLengths& config) {

        capacity_strings_ = capacity_strings - capacity * p + 1; // + 1 to have space for a sentinel after the currently last entry

        strings_ = new char[capacity_strings_];
        string_cursors_ = DacsCapacityArray<dacs_t, builder_t, numSeqs_t>(config.headers_chunk, capacity, config.headers_level,
                                                                          config.headers_extend, config.headers_opt);
        next_string_ = 0;

        prefix_len_ = p;
        prefixes_ = DacsCapacityArray<dacs_t, builder_t, numSeqs_t>(config.prefixes_chunk, capacity, config.prefixes_level,
                                                                    config.prefixes_extend, config.prefixes_opt);

    }

    PrefixDacsStringArray::~PrefixDacsStringArray() {
        delete[] strings_;
    }

    PrefixDacsStringArray::PrefixDacsStringArray(const PrefixDacsStringArray& other) { // copy constructor

        capacity_strings_ = other.capacity_strings_;
        strings_ = new char[capacity_strings_];
        memcpy(strings_, other.strings_, sizeof(char) * capacity_strings_);
        string_cursors_ = other.string_cursors_;
        next_string_ = other.next_string_;

        prefix_len_ = other.prefix_len_;
        pref_num_map_ = other.pref_num_map_;
        num_pref_map_ = other.num_pref_map_;
        prefixes_ = other.prefixes_;

    }

    PrefixDacsStringArray::PrefixDacsStringArray(PrefixDacsStringArray&& other) noexcept { // move constructor

        capacity_strings_ = other.capacity_strings_;
        next_string_ = other.next_string_;
        strings_ = other.strings_; other.strings_ = nullptr;
        string_cursors_ = std::move(other.string_cursors_);

        prefix_len_ = other.prefix_len_;
        pref_num_map_ = other.pref_num_map_;
        num_pref_map_ = other.num_pref_map_;
        prefixes_ = std::move(other.prefixes_);

    }

    PrefixDacsStringArray& PrefixDacsStringArray::operator=(const PrefixDacsStringArray& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] strings_;

        // copy new resources
        capacity_strings_ = other.capacity_strings_;
        strings_ = new char[capacity_strings_];
        memcpy(strings_, other.strings_, sizeof(char) * capacity_strings_);
        string_cursors_ = other.string_cursors_;
        next_string_ = other.next_string_;

        prefix_len_ = other.prefix_len_;
        pref_num_map_ = other.pref_num_map_;
        num_pref_map_ = other.num_pref_map_;
        prefixes_ = other.prefixes_;

        return *this;

    }

    PrefixDacsStringArray& PrefixDacsStringArray::operator=(PrefixDacsStringArray&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] strings_;

        // copy / transfer new resources
        capacity_strings_ = other.capacity_strings_;
        next_string_ = other.next_string_;
        strings_ = other.strings_; other.strings_ = nullptr;
        string_cursors_ = std::move(other.string_cursors_);

        prefix_len_ = other.prefix_len_;
        pref_num_map_ = other.pref_num_map_;
        num_pref_map_ = other.num_pref_map_;
        prefixes_ = std::move(other.prefixes_);

        return *this;

    }

    PrefixDacsStringArray* PrefixDacsStringArray::clone() const { // deep-copy clone method
        return new PrefixDacsStringArray(*this);
    }

    void PrefixDacsStringArray::push_back(const std::string& s) {

        extend(s);

        string_cursors_.push_back(next_string_);
        strcpy(strings_ + next_string_, s.c_str() + prefix_len_);
        next_string_ += s.length() - prefix_len_ + 1;

        std::string prefix = s.substr(0, prefix_len_);
        auto iter = pref_num_map_.find(prefix);

        if (iter == pref_num_map_.end()) { // previously unseen prefix -> extend mapping

            pref_num_map_[prefix] = num_pref_map_.size();
            prefixes_.push_back(num_pref_map_.size());
            num_pref_map_.push_back(prefix);

        } else { // previously seen prefix -> reuse mapping
            prefixes_.push_back(iter->second);
        }

    }

    std::unique_ptr<SequenceWrapper> PrefixDacsStringArray::get(const numSeqs_t i, const lenSeqs_t len) const {
        return std::unique_ptr<SequenceWrapper>(
                new PrefixSeqWrapper(this, i, num_pref_map_[prefixes_.get(i)].length() + strlen(strings_ + string_cursors_.get(i))));
    }

    std::string PrefixDacsStringArray::get_str(const numSeqs_t i, const lenSeqs_t len) const {
        return num_pref_map_[prefixes_.get(i)] + std::string(strings_ + string_cursors_.get(i));
    }

    void PrefixDacsStringArray::extend(const std::string& s) {

        size_t l = s.length();
        unsigned long long len = next_string_;
        if (len + l - prefix_len_ + 1 > capacity_strings_) {

            char* tmp1 = new char[std::max(len + l - prefix_len_ + 1, 2 * capacity_strings_)];
            memcpy(tmp1, strings_, len * sizeof(char));
            delete[] strings_;
            strings_ = tmp1;

            capacity_strings_ = std::max(len + l - prefix_len_ + 1, 2 * capacity_strings_);

        }

    }

    void PrefixDacsStringArray::finalise() {

        string_cursors_.finalise();
        prefixes_.finalise();

    }

    void PrefixDacsStringArray::sort(const std::vector<numSeqs_t>& permutation) {

        string_cursors_.sort(permutation);
        prefixes_.sort(permutation);

    }

    size_t PrefixDacsStringArray::size_in_bytes() const {
        return sizeof(PrefixDacsStringArray) // PrefixDacsStringArray itself
            + sizeof(char) * capacity_strings_ // strings_
            + string_cursors_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) // string_cursors_ (do not count size of DacsCapacityArray object twice)
            + prefixes_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) // prefixes_ (do not count size of DacsCapacityArray object twice)
            + map_str_uchar_size_in_bytes(pref_num_map_) - sizeof(std::map<std::string, unsigned char>) // pref_num_map_ (do not count size of map object twice)
            + vec_str_size_in_bytes(num_pref_map_) - sizeof(std::vector<std::string>); // num_pref_map_ (do not count size of vector object twice)
    }

    void PrefixDacsStringArray::show_memory() const {

        std::cout << "#  * PrefixDacsStringArray " << std::endl;
        std::cout << "#  * " << std::endl;
        std::cout << "#  * sizeof(PrefixDacsStringArray): " << sizeof(PrefixDacsStringArray) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char): " << sizeof(unsigned char) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char*): " << sizeof(unsigned char*) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned long long): " << sizeof(unsigned long long) << " bytes" << std::endl;
        std::cout << "#  * sizeof(unsigned char): " << sizeof(unsigned char) << " bytes" << std::endl;
        std::cout << "#  * sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>): " << sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) << " bytes" << std::endl;
        std::cout << "#  * sizeof(std::map<std::string, unsigned char>): " << sizeof(std::map<std::string, unsigned char>) << " bytes" << std::endl;
        std::cout << "#  * sizeof(std::vector<std::string>): " << sizeof(std::vector<std::string>) << " bytes" << std::endl;
        std::cout << "#  * " << std::endl;
        std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
        std::cout << "#  * strings_: " << (sizeof(char) * capacity_strings_) << " bytes" << std::endl;
        std::cout << "#  * string_cursors_: " << string_cursors_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "#  * prefixes_: " << prefixes_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "#  * pref_num_map_: " << map_str_uchar_size_in_bytes(pref_num_map_) << " bytes" << std::endl;
        std::cout << "#  * num_pref_map_: " << vec_str_size_in_bytes(num_pref_map_) << " bytes" << std::endl;

    }

    /* === PrefixDacsStringArray::PrefixSeqWrapper === */

    PrefixDacsStringArray::PrefixSeqWrapper::PrefixSeqWrapper(const PrefixDacsStringArray* arr, size_t i, lenSeqs_t len)
        : arr_(arr), index_(i), len_(len) {
        // nothing else to do
    }

    char PrefixDacsStringArray::PrefixSeqWrapper::operator[](size_t i) const {
        return arr_->get_char(index_, i);
    }

    lenSeqs_t PrefixDacsStringArray::PrefixSeqWrapper::length() const {
        return len_;
    }

    std::string PrefixDacsStringArray::PrefixSeqWrapper::to_string() const {
        return arr_->get_prefix(index_) + arr_->get_suffix(index_);
    }

    void PrefixDacsStringArray::PrefixSeqWrapper::write(std::ostream& os) const {
        os << arr_->get_prefix(index_) << arr_->get_suffix(index_);
    }

    PrefixDacsStringArray::PrefixSeqWrapper* PrefixDacsStringArray::PrefixSeqWrapper::clone() const {
        return new PrefixSeqWrapper(*this);
    }

    std::unique_ptr<SubstringWrapper> PrefixDacsStringArray::PrefixSeqWrapper::substr(lenSeqs_t start, lenSeqs_t len) const {
        return std::unique_ptr<SubstringWrapper>(new PrefixSubWrapper(arr_, index_, start, len));
    }

    size_t PrefixDacsStringArray::PrefixSeqWrapper::size_in_bytes() const {
        return sizeof(PrefixSeqWrapper);
    }

    bool PrefixDacsStringArray::PrefixSeqWrapper::less(const SequenceWrapper& other) const {

        auto& o = static_cast<const PrefixSeqWrapper&>(other);

        return (arr_->get_prefix(index_) < o.arr_->get_prefix(o.index_))
                || ((arr_->prefixes_.get(index_) == o.arr_->prefixes_.get(o.index_))
                    && (strcmp(arr_->strings_ + arr_->string_cursors_.get(index_), o.arr_->strings_ + o.arr_->string_cursors_.get(o.index_)) < 0));

    }

    bool PrefixDacsStringArray::PrefixSeqWrapper::equals(const SequenceWrapper& other) const {

        auto& o = static_cast<const PrefixSeqWrapper&>(other);

        return (arr_->prefixes_.get(index_) == o.arr_->prefixes_.get(o.index_)) &&
            (strcmp(arr_->strings_ + arr_->string_cursors_.get(index_), o.arr_->strings_ + o.arr_->string_cursors_.get(o.index_)) == 0);

    }

    /* === PrefixDacsStringArray::PrefixSubWrapper === */

    PrefixDacsStringArray::PrefixSubWrapper::PrefixSubWrapper(const PrefixDacsStringArray* arr, size_t i, size_t start, lenSeqs_t len)
        : arr_(arr), index_(i), start_(start), len_(len) {
        // nothing else to do
    }

    char PrefixDacsStringArray::PrefixSubWrapper::operator[](size_t i) const {
        return arr_->get_char(index_, start_ + i);
    }

    lenSeqs_t PrefixDacsStringArray::PrefixSubWrapper::length() const {
        return len_;
    }

    std::string PrefixDacsStringArray::PrefixSubWrapper::to_string() const {
        return arr_->get_string(index_, start_, len_);
    }

    PrefixDacsStringArray::PrefixSubWrapper* PrefixDacsStringArray::PrefixSubWrapper::clone() const {
        return new PrefixSubWrapper(*this);
    }

    void PrefixDacsStringArray::PrefixSubWrapper::shift() {
        start_++;
    }

    size_t PrefixDacsStringArray::PrefixSubWrapper::size_in_bytes() const {
        return sizeof(PrefixSubWrapper);
    }

    bool PrefixDacsStringArray::PrefixSubWrapper::equals(const SubstringWrapper& other) const {

        auto& o = static_cast<const PrefixSubWrapper&>(other);

        bool eq = true;
        for (auto i = 0; eq && i < len_; i++) {
            eq = arr_->get_char(index_, start_ + i) == o[i];
        }

        return eq;

    }




    char rev_acgtu_map[4] = {'A', 'C', 'G', 'T'};

}
