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

#ifndef GEFAST_STRINGARRAYS_HPP
#define GEFAST_STRINGARRAYS_HPP

#include <queue>
#include <sstream>

#include "../Base.hpp"
#include "../modes/SpaceLevenshteinMode.hpp"
#include "Basics.hpp"
#include "NumericArrays.hpp"
#include "dacs/Utility.hpp"

namespace GeFaST {

    /*
     * Basic implementation of StringComponent storing the strings as a C-style array,
     * consecutively containing all strings (including the \0 symbols).
     * Access to the individual strings is provided through a second array of pointers into the string-containing array.
     */
    class CollectiveStringArray : public StringComponent {
    public:
        CollectiveStringArray() = default;

        CollectiveStringArray(const unsigned long long capacity, const unsigned long long capacity_strings);

        ~CollectiveStringArray();

        CollectiveStringArray(const CollectiveStringArray& other);

        CollectiveStringArray(CollectiveStringArray&& other) noexcept;

        CollectiveStringArray& operator=(const CollectiveStringArray& other);

        CollectiveStringArray& operator=(CollectiveStringArray&& other) noexcept;

        CollectiveStringArray* clone() const override;

        void push_back(const std::string& s) override;

        std::unique_ptr<SequenceWrapper> get(const numSeqs_t i, const lenSeqs_t len) const override;

        std::string get_str(const numSeqs_t i, const lenSeqs_t len) const override;

        void extend(const std::string& s) override;

        void finalise() override;

		void sort(const std::vector<numSeqs_t>& permutation) override;

        size_t size_in_bytes() const override;

        void show_memory() const override;

    private:
        char* strings_ = nullptr; // overall string array, each string ends with a \0
        char* next_string_ = nullptr; // position at which the next string would be inserted
        char** string_pointers_ = nullptr; // array of pointers into strings_ marking the start positions of the individual strings
        unsigned long long capacity_strings_ = 0; // capacity of strings_

        numSeqs_t size_ = 0; // current number of amplicons in the collection
        numSeqs_t capacity_ = 0; // capacity of the collection (in terms of the number of amplicons)

    };



    /*
     * Implementation of StringComponent designed for storing headers / identifiers
     * based on the observation that the number of different header prefixes (up to a certain length) is quite small.
     * Therefore, these prefixes are replaced by a single number and only the remaining suffixes are stored
     * as a consecutive C-style array (as in CollectiveStringArray).
     * Access to the individual suffixes is provided through a second array of pointers into the string-containing array.
     */
    class PrefixStringArray : public StringComponent {
    public:
        PrefixStringArray() = default;

        PrefixStringArray(const unsigned char p, const unsigned long long capacity, const unsigned long long capacity_strings);

        ~PrefixStringArray();

        PrefixStringArray(const PrefixStringArray& other);

        PrefixStringArray(PrefixStringArray&& other) noexcept;

        PrefixStringArray& operator=(const PrefixStringArray& other);

        PrefixStringArray& operator=(PrefixStringArray&& other) noexcept;

        PrefixStringArray* clone() const override;

        void push_back(const std::string& s) override;

        std::unique_ptr<SequenceWrapper> get(const numSeqs_t i, const lenSeqs_t len) const override;

        std::string get_str(const numSeqs_t i, const lenSeqs_t len) const override;

        void extend(const std::string& s) override;

        void finalise() override;

        void sort(const std::vector<numSeqs_t>& permutation) override;

        size_t size_in_bytes() const override;

        void show_memory() const override;

    private:
        char* strings_ = nullptr; // overall string (suffix) array, each string ends with a \0
        char* next_string_ = nullptr; // position at which the next string would be inserted
        char** string_pointers_ = nullptr; // array of pointers into strings_ marking the start positions of the individual strings
        unsigned long long capacity_strings_ = 0; // capacity of strings_

        numSeqs_t size_ = 0; // current number of amplicons in the collection
        numSeqs_t capacity_ = 0; // capacity of the collection (in terms of the number of amplicons)

        unsigned char prefix_len_ = 0; // length of the identifier prefixes turned into numbers
        std::map<std::string, unsigned char> pref_num_map_; // mapping between identifier prefix and numeric replacement
        std::vector<std::string> num_pref_map_; // inverse mapping between numeric replacement and identifier prefix
        unsigned char* prefixes_ = nullptr; // array storing the prefix number of each string

        /*
         * Helper methods encapsulating the case distinctions when accessing a single character
         * and extracting a substring, respectively.
         */
        inline char get_char(const size_t index, const size_t pos) const {
            return (pos < prefix_len_) ? num_pref_map_[prefixes_[index]][pos] : string_pointers_[index][pos - prefix_len_];
        }
        inline std::string get_string(const size_t index, const size_t start, const lenSeqs_t len) const {

            std::string tmp;
            if (start < prefix_len_) tmp.append(num_pref_map_[prefixes_[index]].substr(start, len)); // substr() does not extend past the end

            if (tmp.length() != len) { // (remaining) requested characters in the suffix part
                tmp.append(string_pointers_[index] + (start > prefix_len_) * (start - prefix_len_), len - tmp.length());
            }

            return tmp;

        }

        /*
         * Convenience methods for extracting the prefix and suffix of a string.
         */
        inline std::string get_prefix(const size_t index) const {
            return num_pref_map_[prefixes_[index]];
        }
        inline std::string get_suffix(const size_t index) const {
            return string_pointers_[index];
        }


        /*
         * Wrapper for access to whole stored strings.
         */
        class PrefixSeqWrapper : public SequenceWrapper {

        public:
            PrefixSeqWrapper(const PrefixStringArray* arr, size_t i, lenSeqs_t len);

            ~PrefixSeqWrapper() override = default;

            char operator[](size_t i) const override;

            lenSeqs_t length() const override;

            std::string to_string() const override;

            void write(std::ostream& os) const override;

            PrefixSeqWrapper* clone() const override;

            std::unique_ptr<SubstringWrapper> substr(lenSeqs_t start, lenSeqs_t len) const override;

            size_t size_in_bytes() const override;

        private:
            // empty strings are never compared
            bool less(const SequenceWrapper& other) const override;

            bool equals(const SequenceWrapper& other) const override;


            const PrefixStringArray* arr_; // container storing the sequence
            size_t index_; // internal index of the sequence
            lenSeqs_t len_; // length of the sequence

        };

        /*
         * Wrapper for access to substrings.
         */
        class PrefixSubWrapper : public SubstringWrapper {

        public:
            PrefixSubWrapper(const PrefixStringArray* arr, size_t i, size_t start, lenSeqs_t len);

            ~PrefixSubWrapper() override = default;

            char operator[](size_t i) const override;

            lenSeqs_t length() const override;

            std::string to_string() const override;

            PrefixSubWrapper* clone() const override;

            void shift() override;

            size_t size_in_bytes() const override;

        private:
            // only substrings / segments of equal length are compared
            bool equals(const SubstringWrapper& other) const override;

            const PrefixStringArray* arr_; // container storing the sequence
            size_t index_; // internal index of the sequence
            size_t start_; // start position in the sequence
            lenSeqs_t len_; // length of the substring

        };

    };

    /*
     * Variant of PrefixStringArray: the cursors pointing into the overall strings array
     * and the prefix numbers are stored in as DACs sequences.
     */
    class PrefixDacsStringArray : public StringComponent {
    public:
        PrefixDacsStringArray() = default;

        PrefixDacsStringArray(const unsigned char p, const unsigned long long capacity,
                              const unsigned long long capacity_strings, const ChunkLengths& config);

        ~PrefixDacsStringArray();

        PrefixDacsStringArray(const PrefixDacsStringArray& other);

        PrefixDacsStringArray(PrefixDacsStringArray&& other) noexcept;

        PrefixDacsStringArray& operator=(const PrefixDacsStringArray& other);

        PrefixDacsStringArray& operator=(PrefixDacsStringArray&& other) noexcept;

        PrefixDacsStringArray* clone() const override;

        void push_back(const std::string& s) override;

        std::unique_ptr<SequenceWrapper> get(const numSeqs_t i, const lenSeqs_t len) const override;

        std::string get_str(const numSeqs_t i, const lenSeqs_t len) const override;

        void extend(const std::string& s) override;

        void finalise() override;

        void sort(const std::vector<numSeqs_t>& permutation) override;

        size_t size_in_bytes() const override;

        void show_memory() const override;

    private:
        typedef SplitDacs<numSeqs_t, size_t> dacs_t;
        typedef SplitDacsDirectBuilder<numSeqs_t, size_t> builder_t;

        char* strings_ = nullptr; // overall string (suffix) array, each string ends with a \0
        size_t next_string_ = 0; // position at which the next string would be inserted
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> string_cursors_; // array of cursors into strings_ marking the start positions of the individual sequences
        unsigned long long capacity_strings_ = 0; // capacity of strings_

        unsigned char prefix_len_ = 0; // length of the identifier prefixes turned into numbers
        std::map<std::string, unsigned char> pref_num_map_; // mapping between identifier prefix and numeric replacement
        std::vector<std::string> num_pref_map_; // inverse mapping between numeric replacement and identifier prefix
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> prefixes_; // array storing the prefix number of each string

        /*
         * Helper methods encapsulating the case distinctions when accessing a single character
         * and extracting a substring, respectively.
         */
        inline char get_char(const size_t index, const size_t pos) const {
            return (pos < prefix_len_) ? num_pref_map_[prefixes_.get(index)][pos] : (strings_ + string_cursors_.get(index))[pos - prefix_len_];
        }
        inline std::string get_string(const size_t index, const size_t start, const lenSeqs_t len) const {

            std::string tmp;
            if (start < prefix_len_) tmp.append(num_pref_map_[prefixes_.get(index)].substr(start, len)); // substr() does not extend past the end

            if (tmp.length() != len) { // (remaining) requested characters in the suffix part
                tmp.append((strings_ + string_cursors_.get(index)) + (start > prefix_len_) * (start - prefix_len_), len - tmp.length());
            }

            return tmp;

        }

        /*
         * Convenience methods for extracting the prefix and suffix of a string.
         */
        inline std::string get_prefix(const size_t index) const {
            return num_pref_map_[prefixes_.get(index)];
        }
        inline std::string get_suffix(const size_t index) const {
            return strings_ + string_cursors_.get(index);
        }


        /*
         * Wrapper for access to whole stored strings.
         */
        class PrefixSeqWrapper : public SequenceWrapper {

        public:
            PrefixSeqWrapper(const PrefixDacsStringArray* arr, size_t i, lenSeqs_t len);

            ~PrefixSeqWrapper() override = default;

            char operator[](size_t i) const override;

            lenSeqs_t length() const override;

            std::string to_string() const override;

            void write(std::ostream& os) const override;

            PrefixSeqWrapper* clone() const override;

            std::unique_ptr<SubstringWrapper> substr(lenSeqs_t start, lenSeqs_t len) const override;

            size_t size_in_bytes() const override;

        private:
            // empty strings are never compared
            bool less(const SequenceWrapper& other) const override;

            bool equals(const SequenceWrapper& other) const override;


            const PrefixDacsStringArray* arr_; // container storing the sequence
            size_t index_; // internal index of the sequence
            lenSeqs_t len_; // length of the sequence

        };

        /*
         * Wrapper for access to substrings.
         */
        class PrefixSubWrapper : public SubstringWrapper {

        public:
            PrefixSubWrapper(const PrefixDacsStringArray* arr, size_t i, size_t start, lenSeqs_t len);

            ~PrefixSubWrapper() override = default;

            char operator[](size_t i) const override;

            lenSeqs_t length() const override;

            std::string to_string() const override;

            PrefixSubWrapper* clone() const override;

            void shift() override;

            size_t size_in_bytes() const override;

        private:
            // only substrings / segments of equal length are compared
            bool equals(const SubstringWrapper& other) const override;

            const PrefixDacsStringArray* arr_; // container storing the sequence
            size_t index_; // internal index of the sequence
            size_t start_; // start position in the sequence
            lenSeqs_t len_; // length of the substring

        };

    };




    /*
     * Implementation of StringComponent designed for storing headers / identifiers
     * based on the observation that headers are only accessed as a whole or compared from left to right.
     * Before storing the actual headers, a Huffman code is determined for the (previously recorded) distribution
     * of the characters.
     * The encoded strings are stored consecutively (including the encoded \0 symbols) in a bit vector.
     * Access to the individual strings is provided through a second array of cursors pointing at the start positions.
     */
    template<typename S = size_t>
    class HuffmanStringArray : public StringComponent {
    public:
        HuffmanStringArray() = default;

        HuffmanStringArray(const unsigned long long capacity, const std::array<size_t, 256>& char_counts) {

            size_ = 0;
            capacity_ = capacity + 1; // + 1 to have space for a sentinel after the currently last entry

            capacity_strings_ = 1; // + 1 to have space for a sentinel after the currently last entry
            build_huffman_code(char_counts); // build code tree / code words and add require space to capacity_strings_
            capacity_strings_ = (capacity_strings_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); // turn into allocation units

            strings_ = new S[capacity_strings_];
            for (auto i = 0; i < capacity_strings_; i++) {
                strings_[i] = 0;
            }
            string_cursors_ = new size_t[capacity_];
            next_string_ = 0;

        }

        ~HuffmanStringArray() {

            delete[] string_cursors_;
            delete[] strings_;
            delete code_tree_;

        }

        HuffmanStringArray(const HuffmanStringArray& other) { // copy constructor

            capacity_ = other.capacity_;
            size_ = other.size_;

            capacity_strings_ = other.capacity_strings_;
            strings_ = new S[capacity_strings_];
            memcpy(strings_, other.strings_, sizeof(char) * capacity_strings_);
            string_cursors_ = new size_t[capacity_];
            for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
                string_cursors_[i] = other.string_cursors_[i];
            }
            next_string_ = other.next_string_;

            code_tree_ = other.code_tree_->clone();
            code_words_ = other.code_words_;
            code_lengths_ = other.code_lengths_;

        }

        HuffmanStringArray(HuffmanStringArray&& other) noexcept { // move constructor

            capacity_ = other.capacity_;
            size_ = other.size_;

            capacity_strings_ = other.capacity_strings_;
            next_string_ = other.next_string_;
            strings_ = other.strings_; other.strings_ = nullptr;
            string_cursors_ = other.string_cursors_; other.string_cursors_ = nullptr;

            code_tree_ = other.code_tree_; other.code_tree_ = nullptr;
            code_words_ = other.code_words_;
            code_lengths_ = other.code_lengths_;

        }

        HuffmanStringArray& operator=(const HuffmanStringArray& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] string_cursors_;
            delete[] strings_;
            delete code_tree_;

            // copy new resources
            capacity_ = other.capacity_;
            size_ = other.size_;

            capacity_strings_ = other.capacity_strings_;
            strings_ = new S[capacity_strings_];
            memcpy(strings_, other.strings_, sizeof(char) * capacity_strings_);
            string_cursors_ = new size_t[capacity_];
            for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
                string_cursors_[i] = other.string_cursors_[i];
            }
            next_string_ = other.next_string_;

            code_tree_ = other.code_tree_->clone();
            code_words_ = other.code_words_;
            code_lengths_ = other.code_lengths_;

            return *this;

        }

        HuffmanStringArray& operator=(HuffmanStringArray&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] string_cursors_;
            delete[] strings_;
            delete code_tree_;

            // copy / transfer new resources
            capacity_ = other.capacity_;
            size_ = other.size_;

            capacity_strings_ = other.capacity_strings_;
            next_string_ = other.next_string_;
            strings_ = other.strings_; other.strings_ = nullptr;
            string_cursors_ = other.string_cursors_; other.string_cursors_ = nullptr;

            code_tree_ = other.code_tree_; other.code_tree_ = nullptr;
            code_words_ = other.code_words_;
            code_lengths_ = other.code_lengths_;

            return *this;

        }

        HuffmanStringArray* clone() const override { // deep-copy clone method
            return new HuffmanStringArray(*this);
        }

        void push_back(const std::string& s) override {

            extend(s);

            string_cursors_[size_] = next_string_;
            for (auto c : s) {

                write_bits<unsigned char, S>(strings_, next_string_, code_lengths_[c], code_words_[c]);
                next_string_ += code_lengths_[c];

            }
            write_bits<unsigned char, S>(strings_, next_string_, code_lengths_['\0'], code_words_['\0']);
            next_string_ += code_lengths_['\0'];

            size_++;
            string_cursors_[size_] = next_string_; // sentinel for length calculations

        }

        std::unique_ptr<SequenceWrapper> get(const numSeqs_t i, const lenSeqs_t len) const override {
            return std::unique_ptr<SequenceWrapper>(new HuffmanSeqWrapper(this, i));
        }

        std::string get_str(const numSeqs_t i, const lenSeqs_t len) const override {
            return get_string(i);
        }

        void extend(const std::string& s) override {

            size_t len_s = code_lengths_['\0'];
            for (auto c : s) len_s += code_lengths_[c];

            if ((next_string_ + len_s + 1) > capacity_strings_ * (sizeof(S) * 8)) {

                size_t new_capa = (std::max(next_string_ + len_s + 1, 2 * capacity_strings_ * (sizeof(S) * 8)) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8);
                S* tmp1 = new S[new_capa];
                memcpy(tmp1, strings_, capacity_strings_ * sizeof(S));
                delete[] strings_;
                strings_ = tmp1;

                capacity_strings_ = new_capa;

            }

            if (size_ + 1 == capacity_) {

                size_t* tmp1 = new size_t[2 * capacity_];
                memcpy(tmp1, string_cursors_, capacity_ * sizeof(size_t));
                delete[] string_cursors_;
                string_cursors_ = tmp1;

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

                        // swap_values(prev_j, j);
                        std::swap(string_cursors_[prev_j], string_cursors_[j]);

                        done[j] = true;
                        prev_j = j;
                        j = permutation[j];

                    }

                }

            }

        }

        size_t size_in_bytes() const override {
            return sizeof(HuffmanStringArray<S>) // HuffmanStringArray itself
                + sizeof(S) * capacity_strings_ // strings_
                + sizeof(size_t) * capacity_ // string_cursors_
                + tree_size_in_bytes(code_tree_) // code_tree_
                + sizeof(std::array<unsigned char, 256>) // code_words_
                + sizeof(std::array<unsigned char, 256>); // code_lengths_
        }

        void show_memory() const override {

            std::cout << "#  * HuffmanStringArray " << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * sizeof(HuffmanStringArray<S>): " << sizeof(HuffmanStringArray<S>) << " bytes" << std::endl;
            std::cout << "#  * sizeof(S): " << sizeof(S) << " bytes" << std::endl;
            std::cout << "#  * sizeof(S*): " << sizeof(S*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(size_t): " << sizeof(size_t) << " bytes" << std::endl;
            std::cout << "#  * sizeof(size_t*): " << sizeof(size_t*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(numSeqs_t): " << sizeof(numSeqs_t) << " bytes" << std::endl;
            std::cout << "#  * sizeof(HuffNode): " << sizeof(HuffNode) << " bytes" << std::endl;
            std::cout << "#  * sizeof(InternalNode): " << sizeof(InternalNode) << " bytes" << std::endl;
            std::cout << "#  * sizeof(LeafNode): " << sizeof(LeafNode) << " bytes" << std::endl;
            std::cout << "#  * sizeof(std::array<unsigned char, 256>): " << sizeof(std::array<unsigned char, 256>) << " bytes" << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
            std::cout << "#  * strings_: " << (sizeof(S) * capacity_strings_) << " bytes" << std::endl;
            std::cout << "#  * string_cursors_: " << (sizeof(size_t) * capacity_) << " bytes" << std::endl;
            std::cout << "#  * code_tree_: " << tree_size_in_bytes(code_tree_) << " bytes" << std::endl;
            std::cout << "#  * code_words_: " << sizeof(std::array<unsigned char, 256>) << " bytes" << std::endl;
            std::cout << "#  * code_lengths_: " << sizeof(std::array<unsigned char, 256>) << " bytes" << std::endl;

        }

    private:
        /*
         * Base class for nodes in Huffman tree.
         * The priority of the node corresponds the number of occurrences of the "character" represented by the node.
         */
        struct HuffNode {

            size_t priority;

            HuffNode(size_t p) : priority(p) {
                // nothing else to do
            }

            virtual ~HuffNode() = default;

            HuffNode(const HuffNode& other) = default; // copy constructor

            HuffNode(HuffNode&& other) noexcept = default; // move constructor

            HuffNode& operator=(const HuffNode& other) = default; // copy assignment operator

            HuffNode& operator=(HuffNode&& other) noexcept = default; // move assignment operator

            virtual HuffNode* clone() const = 0; // deep-copy clone method

            virtual bool is_leaf() const = 0;

            virtual char get_char() const = 0;

            virtual HuffNode* get_child(unsigned char i) const = 0;

            virtual size_t size_in_bytes() const = 0;

        };

        S* strings_ = nullptr; // overall string array, strings end with a \0-equivalent
        size_t next_string_ = 0; // position at which the next string would be inserted (bit position in strings_)
        size_t* string_cursors_ = nullptr; // array of cursors into strings_ marking the start positions of the individual strings (bit positions)
        size_t capacity_strings_ = 0; // capacity of strings_ (in number of allocated S)

        numSeqs_t size_ = 0; // current number of amplicons in the collection
        numSeqs_t capacity_ = 0; // capacity of the collection (in terms of the number of amplicons)

        HuffNode* code_tree_ = nullptr; // Huffman tree for decoding
        std::array<unsigned char, 256> code_words_; // mapping for encoding characters (value of the code word)
        std::array<unsigned char, 256> code_lengths_; // helper mapping for encoding characters (length of the code word)

        /*
         * Helper methods encapsulating the decoding of a single character and a substring, respectively.
         */
        inline char get_char(const size_t index, const size_t pos) const {

            size_t char_pos = string_cursors_[index];
            for (auto i = 0; i < pos; i++) {
                decode(char_pos);
            }
            return decode(char_pos);

        }
        inline std::string get_string(const size_t index, const size_t start, const lenSeqs_t len) const {

            std::string tmp(len, ' ');
            size_t char_pos = string_cursors_[index];
            for (auto i = 0; i < start; i++) {
                decode(char_pos);
            }
            for (auto i = 0; i < len; i++) {
                tmp[i] = decode(char_pos);
            }

            return tmp;

        }

        /*
         * Helper method decoding a full string.
         */
        inline std::string get_string(const size_t index) const {

            std::string tmp;
            size_t char_pos = string_cursors_[index];
            for (char c = decode(char_pos); c != '\0'; c = decode(char_pos)) {
                tmp.push_back(c);
            }

            return tmp;

        }

        /*
         * Helper method determining the length of a string.
         */
        inline lenSeqs_t get_length(const size_t index) const {

            lenSeqs_t len = 0;
            size_t char_pos = string_cursors_[index];
            for (char c = decode(char_pos); c != '\0'; c = decode(char_pos)) {
                len++;
            }

            return len;

        }

        /*
         * Helper method for decoding a single character starting at the provided position in the bit vector
         * containing all encoded strings.
         */
        char decode(size_t& pos) const {

            HuffNode* node = code_tree_;
            while (!node->is_leaf()) {

                node = node->get_child(get_bit<S>(strings_, pos));
                pos++;

            }

            return node->get_char();

        }

        /*
         * Helper method for decoding a num-th character after the provided position in the bit vector
         * containing all encoded strings.
         */
        char decode_after(size_t& pos, size_t num) const {

            for (auto i = 0; i < num; i++) {
                decode(pos);
            }

            return decode(pos);

        }

        /*
         * Helper method for decoding a substring of the given length starting at the provided position in the bit vector
         * containing all encoded strings.
         */
        std::string decode(size_t& pos, size_t len) const {

            std::string tmp(len, ' ');
            for (auto i = 0; i < len; i++) {
                tmp[i] = decode(pos);
            }

            return tmp;

        }

        /*
         * Helper method for decoding the suffix of a string starting at the provided position in the bit vector
         * containing all encoded strings.
         */
        std::string decode_to_end(size_t& pos) const {

            std::string tmp;
            for (char c = decode(pos); c != '\0'; c = decode(pos)) {
                tmp.push_back(c);
            }

            return tmp;

        }


        // comparison function for nodes in the Huffman tree (priority queue)
        struct greaterHuffNode {
            bool operator()(const HuffNode* lhs, const HuffNode* rhs) const {
                return lhs->priority > rhs->priority;
            }
        };

        /*
         * Internal node in a Huffman tree representing the combined priority of the subsumed characters.
         */
        struct InternalNode : public HuffNode {

            std::array<HuffNode*, 2> children;

            InternalNode(HuffNode* left, HuffNode* right) : HuffNode(left->priority + right->priority) {

                children[0] = left;
                children[1] = right;

            }

            ~InternalNode() {

                delete children[0];
                delete children[1];

            }

            InternalNode(const InternalNode& other) : HuffNode(other) { // copy constructor

                children[0] = other.children[0]->clone();
                children[1] = other.children[1]->clone();

            }

            InternalNode(InternalNode&& other) noexcept : HuffNode(other) { // move constructor

                children[0] = other.children[0]; other.children[0] = nullptr;
                children[1] = other.children[1]; other.children[1] = nullptr;

            }

            InternalNode& operator=(const InternalNode& other) { // copy assignment operator

                // check for self-assignment
                if (&other == this) {
                    return *this;
                }

                HuffNode::operator=(other);

                // release old resources
                delete children[0];
                delete children[1];

                // copy new resources
                children[0] = other.children[0]->clone();
                children[1] = other.children[1]->clone();

                return *this;

            }

            InternalNode& operator=(InternalNode&& other) noexcept { // move assignment operator

                // check for self-assignment
                if (&other == this) {
                    return *this;
                }

                HuffNode::operator=(other);

                // release old resources
                delete children[0];
                delete children[1];

                // copy / transfer new resources
                children[0] = other.children[0]; other.children[0] = nullptr;
                children[1] = other.children[1]; other.children[1] = nullptr;

                return *this;

            }

            InternalNode* clone() const override { // deep-copy clone method
                return new InternalNode(*this);
            }

            bool is_leaf() const override {
                return false;
            }

            char get_char() const override {
                std::cerr << "WARNING: get_char() should never be called on an internal HuffNode." << std::endl;
                return 0;
            }

            HuffNode* get_child(unsigned char i) const override {
                return children[i];
            }

            size_t size_in_bytes() const override {
                return sizeof(InternalNode);
            }

        };

        /*
         * Leaf node in a Huffman tree representing a single character and its priority.
         */
        struct LeafNode : public HuffNode {

            char character;

            LeafNode(char c, size_t p) : HuffNode(p) {
                character = c;
            }

            LeafNode(const LeafNode& other) : HuffNode(other) { // copy constructor
                character = other.character;
            }

            LeafNode(LeafNode&& other) noexcept : HuffNode(other) { // move constructor
                character = other.character;
            }

            LeafNode& operator=(const LeafNode& other) { // copy assignment operator

                // check for self-assignment
                if (&other == this) {
                    return *this;
                }

                HuffNode::operator=(other);

                character = other.character;

                return *this;

            }

            LeafNode& operator=(LeafNode&& other) noexcept { // move assignment operator

                // check for self-assignment
                if (&other == this) {
                    return *this;
                }

                HuffNode::operator=(other);

                character = other.character;

                return *this;

            }

            LeafNode* clone() const override { // deep-copy clone method
                return new LeafNode(*this);
            }

            bool is_leaf() const override {
                return true;
            }

            char get_char() const override {
                return character;
            }

            HuffNode* get_child(unsigned char i) const override {
                std::cerr << "WARNING: get_child() should never be called on a leaf HuffNode." << std::endl;
                return nullptr;
            }

            size_t size_in_bytes() const override {
                return sizeof(LeafNode);
            }

        };

        /*
         * Build Huffman tree / code by repeatedly combining the nodes with the smallest priority.
         */
        void build_huffman_code(const std::array<size_t, 256>& counts) {

            std::priority_queue<HuffNode*, std::vector<HuffNode*>, greaterHuffNode> pqueue;

            for (auto c = 0; c < 256; c++) {
                if(counts[c] != 0) pqueue.push(new LeafNode(static_cast<char>(c), counts[c]));
            }

            while (pqueue.size() > 1) {

                auto left = pqueue.top(); pqueue.pop();
                auto right = pqueue.top(); pqueue.pop();
                pqueue.push(new InternalNode(left, right));

            }

            code_tree_ = pqueue.top();

            set_code_words(code_tree_, 0, 0);

        }

        /*
         * Determine the code words and their lengths by traversing the Huffman tree.
         */
        void set_code_words(const HuffNode* node, unsigned char prefix, size_t prefix_len) {

            if (node->is_leaf()) {

                code_words_[node->get_char()] = prefix;
                code_lengths_[node->get_char()] = prefix_len;
                capacity_strings_ += node->priority * prefix_len;

            } else {

                set_code_words(node->get_child(0), prefix << 1, prefix_len + 1);
                set_code_words(node->get_child(1), (prefix << 1) | 1, prefix_len + 1);

            }

        }

        /*
         * Determine the memory consumption (in bytes) of the tree with the given node as its root.
         */
        size_t tree_size_in_bytes(HuffNode* root) const {
            return root->size_in_bytes() + ((!root->is_leaf()) ? (tree_size_in_bytes(root->get_child(0)) + tree_size_in_bytes(root->get_child(1))) : 0);
        }


        /*
         * Wrapper for access to whole stored strings.
         */
        class HuffmanSeqWrapper : public SequenceWrapper {

        public:
            HuffmanSeqWrapper(const HuffmanStringArray* arr, size_t i) : arr_(arr), index_(i) {
                // nothing else to do
            }

            ~HuffmanSeqWrapper() override = default;

            char operator[](size_t i) const override {
                return arr_->get_char(index_, i);
            }

            lenSeqs_t length() const override {
                return arr_->get_length(index_);
            }

            std::string to_string() const override {
                return arr_->get_string(index_);
            }

            void write(std::ostream& os) const override {
                os << arr_->get_string(index_);
            }

            HuffmanSeqWrapper* clone() const override {
                return new HuffmanSeqWrapper(*this);
            }

            std::unique_ptr<SubstringWrapper> substr(lenSeqs_t start, lenSeqs_t len) const override {

                size_t char_pos = get_start();
                for (auto i = 0; i < start; i++) {
                    decode(char_pos);
                }

                return std::unique_ptr<SubstringWrapper>(new HuffmanSubWrapper(arr_, char_pos, len));

            }

            size_t size_in_bytes() const override {
                return sizeof(HuffmanSeqWrapper);
            }

            size_t get_start() const {
                return arr_->string_cursors_[index_];
            }

            char decode(size_t& char_pos) const {
                return arr_->decode(char_pos);
            }

        private:
            bool less(const SequenceWrapper& other) const override {

                auto& o = static_cast<const HuffmanSeqWrapper&>(other);

                bool eq = true;
                size_t char_pos = get_start();
                size_t o_char_pos = o.get_start();
                char last_char, o_last_char;
                do {

                    last_char = decode(char_pos);
                    o_last_char = o.decode(o_char_pos);
                    eq = last_char == o_last_char;

                } while (eq && last_char != '\0' && o_last_char != '\0');

                return (eq && last_char < o_last_char) || (last_char < o_last_char);

            }

            bool equals(const SequenceWrapper& other) const override {

                auto& o = static_cast<const HuffmanSeqWrapper&>(other);

                bool eq = true;
                size_t char_pos = get_start();
                size_t o_char_pos = o.get_start();
                char last_char, o_last_char;
                do {

                    last_char = decode(char_pos);
                    o_last_char = o.decode(o_char_pos);
                    eq = last_char == o_last_char;

                } while (eq && last_char != '\0' && o_last_char != '\0');

                return eq && (last_char == o_last_char);

            }


            const HuffmanStringArray* arr_; // container storing the sequence
            size_t index_; // internal index of the sequence

        };

        /*
         * Wrapper for access to substrings.
         */
        class HuffmanSubWrapper : public SubstringWrapper {

        public:
            HuffmanSubWrapper(const HuffmanStringArray* arr, size_t start, lenSeqs_t len) : arr_(arr), start_(start), len_(len) {
                // nothing else to do
            }

            ~HuffmanSubWrapper() override = default;

            char operator[](size_t i) const override {

                size_t pos = start_;
                return arr_->decode_after(pos, i);

            }

            lenSeqs_t length() const override {
                return len_;
            }

            std::string to_string() const override {

                size_t pos = start_;
                return arr_->decode(pos, len_);

            }

            HuffmanSubWrapper* clone() const override {
                return new HuffmanSubWrapper(*this);
            }

            void shift() override {
                arr_->decode(start_);
            }

            size_t size_in_bytes() const override {
                return sizeof(HuffmanSubWrapper);
            }

            size_t get_start() const {
                return start_;
            }

            char decode(size_t& char_pos) const {
                return arr_->decode(char_pos);
            }

        private:
            // only substrings / segments of equal length are compared
            bool equals(const SubstringWrapper& other) const override {

                auto& o = static_cast<const HuffmanSubWrapper&>(other);

                bool eq = true;
                size_t char_pos = get_start();
                size_t o_char_pos = o.get_start();
                for (auto i = 0; eq && i < len_; i++) {
                    eq = decode(char_pos) == o.decode(o_char_pos);
                }

                return eq;

            }

            const HuffmanStringArray* arr_; // container storing the sequence
            size_t start_; // start position in the sequence
            lenSeqs_t len_; // length of the substring

        };

    };

    /*
     * Variant of HuffmanStringArray: the cursors pointing into the overall strings array are stored as a DACs sequence.
     */
    template<typename S = size_t>
    class HuffmanDacsStringArray : public StringComponent {
    public:
        HuffmanDacsStringArray() = default;

        HuffmanDacsStringArray(const unsigned long long capacity, const std::array<size_t, 256>& char_counts, const ChunkLengths& config) {

            capacity_strings_ = 1; // + 1 to have space for a sentinel after the currently last entry
            build_huffman_code(char_counts); // build code tree / code words and add require space to capacity_strings_
            capacity_strings_ = (capacity_strings_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); // turn into allocation units

            strings_ = new S[capacity_strings_];
            for (auto i = 0; i < capacity_strings_; i++) {
                strings_[i] = 0;
            }
            string_cursors_ = DacsCapacityArray<dacs_t, builder_t, numSeqs_t>(config.headers_chunk, capacity, config.headers_level,
                                                                              config.headers_extend, config.headers_opt);
            next_string_ = 0;

        }

        ~HuffmanDacsStringArray() {

            delete[] strings_;
            delete code_tree_;

        }

        HuffmanDacsStringArray(const HuffmanDacsStringArray& other) { // copy constructor

            capacity_strings_ = other.capacity_strings_;
            strings_ = new S[capacity_strings_];
            memcpy(strings_, other.strings_, sizeof(char) * capacity_strings_);
            string_cursors_ = other.string_cursors_;
            next_string_ = other.next_string_;

            code_tree_ = other.code_tree_->clone();
            code_words_ = other.code_words_;
            code_lengths_ = other.code_lengths_;

        }

        HuffmanDacsStringArray(HuffmanDacsStringArray&& other) noexcept { // move constructor

            capacity_strings_ = other.capacity_strings_;
            next_string_ = other.next_string_;
            strings_ = other.strings_; other.strings_ = nullptr;
            string_cursors_ = std::move(other.string_cursors_);

            code_tree_ = other.code_tree_; other.code_tree_ = nullptr;
            code_words_ = other.code_words_;
            code_lengths_ = other.code_lengths_;

        }

        HuffmanDacsStringArray& operator=(const HuffmanDacsStringArray& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] strings_;
            delete code_tree_;

            // copy new resources
            capacity_strings_ = other.capacity_strings_;
            strings_ = new S[capacity_strings_];
            memcpy(strings_, other.strings_, sizeof(char) * capacity_strings_);
            string_cursors_ = other.string_cursors_;
            next_string_ = other.next_string_;

            code_tree_ = other.code_tree_->clone();
            code_words_ = other.code_words_;
            code_lengths_ = other.code_lengths_;

            return *this;

        }

        HuffmanDacsStringArray& operator=(HuffmanDacsStringArray&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] strings_;
            delete code_tree_;

            // copy / transfer new resources
            capacity_strings_ = other.capacity_strings_;
            next_string_ = other.next_string_;
            strings_ = other.strings_; other.strings_ = nullptr;
            string_cursors_ = std::move(other.string_cursors_);

            code_tree_ = other.code_tree_; other.code_tree_ = nullptr;
            code_words_ = other.code_words_;
            code_lengths_ = other.code_lengths_;

            return *this;

        }

        HuffmanDacsStringArray* clone() const override { // deep-copy clone method
            return new HuffmanDacsStringArray(*this);
        }

        void push_back(const std::string& s) override {

            extend(s);

            string_cursors_.push_back(next_string_);
            for (auto c : s) {

                write_bits<unsigned char, S>(strings_, next_string_, code_lengths_[c], code_words_[c]);
                next_string_ += code_lengths_[c];

            }
            write_bits<unsigned char, S>(strings_, next_string_, code_lengths_['\0'], code_words_['\0']);
            next_string_ += code_lengths_['\0'];

        }

        std::unique_ptr<SequenceWrapper> get(const numSeqs_t i, const lenSeqs_t len) const override {
            return std::unique_ptr<SequenceWrapper>(new HuffmanSeqWrapper(this, i));
        }

        std::string get_str(const numSeqs_t i, const lenSeqs_t len) const override {
            return get_string(i);
        }

        void extend(const std::string& s) override {

            size_t len_s = code_lengths_['\0'];
            for (auto c : s) len_s += code_lengths_[c];

            if ((next_string_ + len_s + 1) > capacity_strings_ * (sizeof(S) * 8)) {

                size_t new_capa = (std::max(next_string_ + len_s + 1, 2 * capacity_strings_ * (sizeof(S) * 8)) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8);
                S* tmp1 = new S[new_capa];
                memcpy(tmp1, strings_, capacity_strings_ * sizeof(S));
                delete[] strings_;
                strings_ = tmp1;

                capacity_strings_ = new_capa;

            }

        }

        void finalise() override {
            string_cursors_.finalise();
        }

        void sort(const std::vector<numSeqs_t>& permutation) override {
            string_cursors_.sort(permutation);
        }

        size_t size_in_bytes() const override {
            return sizeof(HuffmanDacsStringArray<S>) // HuffmanDacsStringArray itself
                + sizeof(S) * capacity_strings_ // strings_
                + string_cursors_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) // string_cursors_ (do not count size of DacsCapacityArray object twice)
                + tree_size_in_bytes(code_tree_) // code_tree_
                + sizeof(std::array<unsigned char, 256>) // code_words_
                + sizeof(std::array<unsigned char, 256>); // code_lengths_
        }

        void show_memory() const override {

            std::cout << "#  * HuffmanDacsStringArray " << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * sizeof(HuffmanDacsStringArray<S>): " << sizeof(HuffmanDacsStringArray<S>) << " bytes" << std::endl;
            std::cout << "#  * sizeof(S): " << sizeof(S) << " bytes" << std::endl;
            std::cout << "#  * sizeof(S*): " << sizeof(S*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(size_t): " << sizeof(size_t) << " bytes" << std::endl;
            std::cout << "#  * sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>): " << sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) << " bytes" << std::endl;
            std::cout << "#  * sizeof(numSeqs_t): " << sizeof(numSeqs_t) << " bytes" << std::endl;
            std::cout << "#  * sizeof(HuffNode): " << sizeof(HuffNode) << " bytes" << std::endl;
            std::cout << "#  * sizeof(InternalNode): " << sizeof(InternalNode) << " bytes" << std::endl;
            std::cout << "#  * sizeof(LeafNode): " << sizeof(LeafNode) << " bytes" << std::endl;
            std::cout << "#  * sizeof(std::array<unsigned char, 256>): " << sizeof(std::array<unsigned char, 256>) << " bytes" << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
            std::cout << "#  * strings_: " << (sizeof(S) * capacity_strings_) << " bytes" << std::endl;
            std::cout << "#  * string_cursors_: " << string_cursors_.size_in_bytes() << " bytes" << std::endl;
            std::cout << "#  * code_tree_: " << tree_size_in_bytes(code_tree_) << " bytes" << std::endl;
            std::cout << "#  * code_words_: " << sizeof(std::array<unsigned char, 256>) << " bytes" << std::endl;
            std::cout << "#  * code_lengths_: " << sizeof(std::array<unsigned char, 256>) << " bytes" << std::endl;

        }

    private:
        /*
         * Base class for nodes in Huffman tree.
         * The priority of the node corresponds the number of occurrences of the "character" represented by the node.
         */
        struct HuffNode {

            size_t priority;

            HuffNode(size_t p) : priority(p) {
                // nothing else to do
            }

            virtual ~HuffNode() = default;

            HuffNode(const HuffNode& other) = default; // copy constructor

            HuffNode(HuffNode&& other) noexcept = default; // move constructor

            HuffNode& operator=(const HuffNode& other) = default; // copy assignment operator

            HuffNode& operator=(HuffNode&& other) noexcept = default; // move assignment operator

            virtual HuffNode* clone() const = 0; // deep-copy clone method

            virtual bool is_leaf() const = 0;

            virtual char get_char() const = 0;

            virtual HuffNode* get_child(unsigned char i) const = 0;

            virtual size_t size_in_bytes() const = 0;

        };

        typedef SplitDacs<numSeqs_t, size_t> dacs_t;
        typedef SplitDacsDirectBuilder<numSeqs_t, size_t> builder_t;

        S* strings_ = nullptr; // overall string array, strings end with a \0-equivalent
        size_t next_string_ = 0; // position at which the next string would be inserted (bit position in strings_)
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> string_cursors_; // array of cursors into strings_ marking the start positions of the individual sequences (bit positions)
        size_t capacity_strings_ = 0; // capacity of strings_ (in number of allocated S)

        HuffNode* code_tree_ = nullptr; // Huffman tree for decoding
        std::array<unsigned char, 256> code_words_; // mapping for encoding characters (value of the code word)
        std::array<unsigned char, 256> code_lengths_; // helper mapping for encoding characters (length of the code word)

        /*
         * Helper methods encapsulating the decoding of a single character and a substring, respectively.
         */
        inline char get_char(const size_t index, const size_t pos) const {

            size_t char_pos = string_cursors_.get(index);
            for (auto i = 0; i < pos; i++) {
                decode(char_pos);
            }
            return decode(char_pos);

        }
        inline std::string get_string(const size_t index, const size_t start, const lenSeqs_t len) const {

            std::string tmp(len, ' ');
            size_t char_pos = string_cursors_.get(index);
            for (auto i = 0; i < start; i++) {
                decode(char_pos);
            }
            for (auto i = 0; i < len; i++) {
                tmp[i] = decode(char_pos);
            }

            return tmp;

        }

        /*
         * Helper method decoding a full string.
         */
        inline std::string get_string(const size_t index) const {

            std::string tmp;
            size_t char_pos = string_cursors_.get(index);
            for (char c = decode(char_pos); c != '\0'; c = decode(char_pos)) {
                tmp.push_back(c);
            }

            return tmp;

        }

        /*
         * Helper method determining the length of a string.
         */
        inline lenSeqs_t get_length(const size_t index) const {

            lenSeqs_t len = 0;
            size_t char_pos = string_cursors_.get(index);
            for (char c = decode(char_pos); c != '\0'; c = decode(char_pos)) {
                len++;
            }

            return len;

        }

        /*
         * Helper method for decoding a single character starting at the provided position in the bit vector
         * containing all encoded strings.
         */
        char decode(size_t& pos) const {

            HuffNode* node = code_tree_;
            while (!node->is_leaf()) {

                node = node->get_child(get_bit<S>(strings_, pos));
                pos++;

            }

            return node->get_char();

        }

        /*
         * Helper method for decoding a num-th character after the provided position in the bit vector
         * containing all encoded strings.
         */
        char decode_after(size_t& pos, size_t num) const {

            for (auto i = 0; i < num; i++) {
                decode(pos);
            }

            return decode(pos);

        }

        /*
         * Helper method for decoding a substring of the given length starting at the provided position in the bit vector
         * containing all encoded strings.
         */
        std::string decode(size_t& pos, size_t len) const {

            std::string tmp(len, ' ');
            for (auto i = 0; i < len; i++) {
                tmp[i] = decode(pos);
            }

            return tmp;

        }

        /*
         * Helper method for decoding the suffix of a string starting at the provided position in the bit vector
         * containing all encoded strings.
         */
        std::string decode_to_end(size_t& pos) const {

            std::string tmp;
            for (char c = decode(pos); c != '\0'; c = decode(pos)) {
                tmp.push_back(c);
            }

            return tmp;

        }


        // comparison function for nodes in the Huffman tree (priority queue)
        struct greaterHuffNode {
            bool operator()(const HuffNode* lhs, const HuffNode* rhs) const {
                return lhs->priority > rhs->priority;
            }
        };

        /*
         * Internal node in a Huffman tree representing the combined priority of the subsumed characters.
         */
        struct InternalNode : public HuffNode {

            std::array<HuffNode*, 2> children;

            InternalNode(HuffNode* left, HuffNode* right) : HuffNode(left->priority + right->priority) {

                children[0] = left;
                children[1] = right;

            }

            ~InternalNode() {

                delete children[0];
                delete children[1];

            }

            InternalNode(const InternalNode& other) : HuffNode(other) { // copy constructor

                children[0] = other.children[0]->clone();
                children[1] = other.children[1]->clone();

            }

            InternalNode(InternalNode&& other) noexcept : HuffNode(other) { // move constructor

                children[0] = other.children[0]; other.children[0] = nullptr;
                children[1] = other.children[1]; other.children[1] = nullptr;

            }

            InternalNode& operator=(const InternalNode& other) { // copy assignment operator

                // check for self-assignment
                if (&other == this) {
                    return *this;
                }

                HuffNode::operator=(other);

                // release old resources
                delete children[0];
                delete children[1];

                // copy new resources
                children[0] = other.children[0]->clone();
                children[1] = other.children[1]->clone();

                return *this;

            }

            InternalNode& operator=(InternalNode&& other) noexcept { // move assignment operator

                // check for self-assignment
                if (&other == this) {
                    return *this;
                }

                HuffNode::operator=(other);

                // release old resources
                delete children[0];
                delete children[1];

                // copy / transfer new resources
                children[0] = other.children[0]; other.children[0] = nullptr;
                children[1] = other.children[1]; other.children[1] = nullptr;

                return *this;

            }

            InternalNode* clone() const override { // deep-copy clone method
                return new InternalNode(*this);
            }

            bool is_leaf() const override {
                return false;
            }

            char get_char() const override {
                std::cerr << "WARNING: get_char() should never be called on an internal HuffNode." << std::endl;
                return 0;
            }

            HuffNode* get_child(unsigned char i) const override {
                return children[i];
            }

            size_t size_in_bytes() const override {
                return sizeof(InternalNode);
            }

        };

        /*
         * Leaf node in a Huffman tree representing a single character and its priority.
         */
        struct LeafNode : public HuffNode {

            char character;

            LeafNode(char c, size_t p) : HuffNode(p) {
                character = c;
            }

            LeafNode(const LeafNode& other) : HuffNode(other) { // copy constructor
                character = other.character;
            }

            LeafNode(LeafNode&& other) noexcept : HuffNode(other) { // move constructor
                character = other.character;
            }

            LeafNode& operator=(const LeafNode& other) { // copy assignment operator

                // check for self-assignment
                if (&other == this) {
                    return *this;
                }

                HuffNode::operator=(other);

                character = other.character;

                return *this;

            }

            LeafNode& operator=(LeafNode&& other) noexcept { // move assignment operator

                // check for self-assignment
                if (&other == this) {
                    return *this;
                }

                HuffNode::operator=(other);

                character = other.character;

                return *this;

            }

            LeafNode* clone() const override { // deep-copy clone method
                return new LeafNode(*this);
            }

            bool is_leaf() const override {
                return true;
            }

            char get_char() const override {
                return character;
            }

            HuffNode* get_child(unsigned char i) const override {
                std::cerr << "WARNING: get_child() should never be called on a leaf HuffNode." << std::endl;
                return nullptr;
            }

            size_t size_in_bytes() const override {
                return sizeof(LeafNode);
            }

        };

        /*
         * Build Huffman tree / code by repeatedly combining the nodes with the smallest priority.
         */
        void build_huffman_code(const std::array<size_t, 256>& counts) {

            std::priority_queue<HuffNode*, std::vector<HuffNode*>, greaterHuffNode> pqueue;

            for (auto c = 0; c < 256; c++) {
                if(counts[c] != 0) pqueue.push(new LeafNode(static_cast<char>(c), counts[c]));
            }

            while (pqueue.size() > 1) {

                auto left = pqueue.top(); pqueue.pop();
                auto right = pqueue.top(); pqueue.pop();
                pqueue.push(new InternalNode(left, right));

            }

            code_tree_ = pqueue.top();

            set_code_words(code_tree_, 0, 0);

        }

        /*
         * Determine the code words and their lengths by traversing the Huffman tree.
         */
        void set_code_words(const HuffNode* node, unsigned char prefix, size_t prefix_len) {

            if (node->is_leaf()) {

                code_words_[node->get_char()] = prefix;
                code_lengths_[node->get_char()] = prefix_len;
                capacity_strings_ += node->priority * prefix_len;

            } else {

                set_code_words(node->get_child(0), prefix << 1, prefix_len + 1);
                set_code_words(node->get_child(1), (prefix << 1) | 1, prefix_len + 1);

            }

        }

        /*
         * Determine the memory consumption (in bytes) of the tree with the given node as its root.
         */
        size_t tree_size_in_bytes(HuffNode* root) const {
            return root->size_in_bytes() + ((!root->is_leaf()) ? (tree_size_in_bytes(root->get_child(0)) + tree_size_in_bytes(root->get_child(1))) : 0);
        }


        /*
         * Wrapper for access to whole stored strings.
         */
        class HuffmanSeqWrapper : public SequenceWrapper {

        public:
            HuffmanSeqWrapper(const HuffmanDacsStringArray* arr, size_t i) : arr_(arr), index_(i) {
                // nothing else to do
            }

            ~HuffmanSeqWrapper() override = default;

            char operator[](size_t i) const override {
                return arr_->get_char(index_, i);
            }

            lenSeqs_t length() const override {
                return arr_->get_length(index_);
            }

            std::string to_string() const override {
                return arr_->get_string(index_);
            }

            void write(std::ostream& os) const override {
                os << arr_->get_string(index_);
            }

            HuffmanSeqWrapper* clone() const override {
                return new HuffmanSeqWrapper(*this);
            }

            std::unique_ptr<SubstringWrapper> substr(lenSeqs_t start, lenSeqs_t len) const override {

                size_t char_pos = get_start();
                for (auto i = 0; i < start; i++) {
                    decode(char_pos);
                }

                return std::unique_ptr<SubstringWrapper>(new HuffmanSubWrapper(arr_, char_pos, len));

            }

            size_t size_in_bytes() const override {
                return sizeof(HuffmanSeqWrapper);
            }

            size_t get_start() const {
                return arr_->string_cursors_.get(index_);
            }

            char decode(size_t& char_pos) const {
                return arr_->decode(char_pos);
            }

        private:
            bool less(const SequenceWrapper& other) const override {

                auto& o = static_cast<const HuffmanSeqWrapper&>(other);

                bool eq = true;
                size_t char_pos = get_start();
                size_t o_char_pos = o.get_start();
                char last_char, o_last_char;
                do {

                    last_char = decode(char_pos);
                    o_last_char = o.decode(o_char_pos);
                    eq = last_char == o_last_char;

                } while (eq && last_char != '\0' && o_last_char != '\0');

                return (eq && last_char < o_last_char) || (last_char < o_last_char);

            }

            bool equals(const SequenceWrapper& other) const override {

                auto& o = static_cast<const HuffmanSeqWrapper&>(other);

                bool eq = true;
                size_t char_pos = get_start();
                size_t o_char_pos = o.get_start();
                char last_char, o_last_char;
                do {

                    last_char = decode(char_pos);
                    o_last_char = o.decode(o_char_pos);
                    eq = last_char == o_last_char;

                } while (eq && last_char != '\0' && o_last_char != '\0');

                return eq && (last_char == o_last_char);

            }


            const HuffmanDacsStringArray* arr_; // container storing the sequence
            size_t index_; // internal index of the sequence

        };

        /*
         * Wrapper for access to substrings.
         */
        class HuffmanSubWrapper : public SubstringWrapper {

        public:
            HuffmanSubWrapper(const HuffmanDacsStringArray* arr, size_t start, lenSeqs_t len) : arr_(arr), start_(start), len_(len) {
                // nothing else to do
            }

            ~HuffmanSubWrapper() override = default;

            char operator[](size_t i) const override {

                size_t pos = start_;
                return arr_->decode_after(pos, i);

            }

            lenSeqs_t length() const override {
                return len_;
            }

            std::string to_string() const override {

                size_t pos = start_;
                return arr_->decode(pos, len_);

            }

            HuffmanSubWrapper* clone() const override {
                return new HuffmanSubWrapper(*this);
            }

            void shift() override {
                arr_->decode(start_);
            }

            size_t size_in_bytes() const override {
                return sizeof(HuffmanSubWrapper);
            }

            size_t get_start() const {
                return start_;
            }

            char decode(size_t& char_pos) const {
                return arr_->decode(char_pos);
            }

        private:
            // only substrings / segments of equal length are compared
            bool equals(const SubstringWrapper& other) const override {

                auto& o = static_cast<const HuffmanSubWrapper&>(other);

                bool eq = true;
                size_t char_pos = get_start();
                size_t o_char_pos = o.get_start();
                for (auto i = 0; eq && i < len_; i++) {
                    eq = decode(char_pos) == o.decode(o_char_pos);
                }

                return eq;

            }

            const HuffmanDacsStringArray* arr_; // container storing the sequence
            size_t start_; // start position in the sequence
            lenSeqs_t len_; // length of the substring

        };

    };
	
	
	
	// reverse mapping of integers onto nucleotides (1 -> A, 2 -> C, 3 -> G, 4 -> T)
    // effectively turns u/U into T (same behaviour as in Swarm)
    extern char rev_acgtu_map[4];

    /*
     * Implementation of StringComponent designed for storing nucleotide sequences (alphabet of size 4)
     * using 2 instead of 8 bits per character.
     * The encoded strings are stored consecutively (without the \0 symbols) in a bit vector.
     * Access to the individual strings is provided through a second array of cursors pointing at the start positions.
     */
    template<typename S = size_t>
    class TwoBitStringArray : public StringComponent {
    public:
        TwoBitStringArray() = default;

        TwoBitStringArray(const unsigned long long capacity, const unsigned long long capacity_strings) {

            size_ = 0;
            capacity_ = capacity + 1; // + 1 to have space for a sentinel after the currently last entry
            capacity_strings_ = (2 * (capacity_strings + 1 - capacity) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); // + 1 for sentinel, - capacity because final \0 are not stored

            strings_ = new S[capacity_strings_];
            for (auto i = 0; i < capacity_strings_; i++) {
                strings_[i] = 0;
            }
            string_cursors_ = new size_t[capacity_];
            next_string_ = 0;

        }

        ~TwoBitStringArray() {

            delete[] string_cursors_;
            delete[] strings_;

        }

        TwoBitStringArray(const TwoBitStringArray& other) { // copy constructor

            capacity_ = other.capacity_;
            size_ = other.size_;

            capacity_strings_ = other.capacity_strings_;
            strings_ = new S[capacity_strings_];
            memcpy(strings_, other.strings_, sizeof(S) * capacity_strings_);
            string_cursors_ = new size_t[capacity_];
            for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
                string_cursors_[i] = other.string_cursors_[i];
            }
            next_string_ = other.next_string_;

        }

        TwoBitStringArray(TwoBitStringArray&& other) noexcept { // move constructor

            capacity_ = other.capacity_;
            size_ = other.size_;

            capacity_strings_ = other.capacity_strings_;
            next_string_ = other.next_string_;
            strings_ = other.strings_; other.strings_ = nullptr;
            string_cursors_ = other.string_cursors_; other.string_cursors_ = nullptr;

        }

        TwoBitStringArray& operator=(const TwoBitStringArray& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] string_cursors_;
            delete[] strings_;

            // copy new resources
            capacity_ = other.capacity_;
            size_ = other.size_;

            capacity_strings_ = other.capacity_strings_;
            strings_ = new S[capacity_strings_];
            memcpy(strings_, other.strings_, sizeof(S) * capacity_strings_);
            string_cursors_ = new size_t[capacity_];
            for (auto i = 0; i <= size_; i++) { // <= to copy the sentinel
                string_cursors_[i] = other.string_cursors_[i];
            }
            next_string_ = other.next_string_;

            return *this;

        }

        TwoBitStringArray& operator=(TwoBitStringArray&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] string_cursors_;
            delete[] strings_;

            // copy / transfer new resources
            capacity_ = other.capacity_;
            size_ = other.size_;

            capacity_strings_ = other.capacity_strings_;
            next_string_ = other.next_string_;
            strings_ = other.strings_; other.strings_ = nullptr;
            string_cursors_ = other.string_cursors_; other.string_cursors_ = nullptr;

            return *this;

        }

        TwoBitStringArray* clone() const override { // deep-copy clone method
            return new TwoBitStringArray(*this);
        }

        inline void push_back(const std::string& s) override {

            extend(s);

            string_cursors_[size_] = next_string_;
            for (size_t i = 0, inner = next_string_ % (sizeof(S) * 8), outer = next_string_ / (sizeof(S) * 8);
                 i < s.length();
                 i++, inner = (inner + 2) % (sizeof(S) * 8), outer += (inner == 0)) {

                strings_[outer] |= static_cast<S>(acgtu_map[s[i]] - 1) << ((sizeof(S) * 8) - 2 - inner);

            }
            next_string_ += 2 * s.length();

            size_++;
            string_cursors_[size_] = next_string_; // sentinel for length calculations

        }

        inline std::unique_ptr<SequenceWrapper> get(const numSeqs_t i, const lenSeqs_t len) const override {
            return std::unique_ptr<SequenceWrapper>(new TwoBitSeqWrapper(this, string_cursors_[i], len));
        }

        std::string get_str(const numSeqs_t i, const lenSeqs_t len) const override {
            return get_string(string_cursors_[i], len);
        }

        inline void extend(const std::string& s) override {

            size_t len_s = 2 * s.length();

            if ((next_string_ + len_s + 2) > capacity_strings_ * (sizeof(S) * 8)) {

                size_t new_capa = (std::max(next_string_ + len_s + 2, 2 * capacity_strings_ * (sizeof(S) * 8)) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8);
                S* tmp1 = new S[new_capa];
                memcpy(tmp1, strings_, capacity_strings_ * sizeof(S));
                delete[] strings_;
                strings_ = tmp1;

                capacity_strings_ = new_capa;

            }

            if (size_ + 1 == capacity_) {

                size_t* tmp1 = new size_t[2 * capacity_];
                memcpy(tmp1, string_cursors_, capacity_ * sizeof(size_t));
                delete[] string_cursors_;
                string_cursors_ = tmp1;

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

                        // swap_values(prev_j, j);
                        std::swap(string_cursors_[prev_j], string_cursors_[j]);

                        done[j] = true;
                        prev_j = j;
                        j = permutation[j];

                    }

                }

            }

        }

        size_t size_in_bytes() const override {
            return sizeof(TwoBitStringArray<S>) // TwoBitStringArray itself
                + sizeof(S) * capacity_strings_ // strings_
                + sizeof(size_t) * capacity_; // string_cursors_
        }

        void show_memory() const override {

            std::cout << "#  * TwoBitStringArray " << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * sizeof(TwoBitStringArray<S>): " << sizeof(TwoBitStringArray<S>) << " bytes" << std::endl;
            std::cout << "#  * sizeof(S): " << sizeof(S) << " bytes" << std::endl;
            std::cout << "#  * sizeof(S*): " << sizeof(S*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(size_t): " << sizeof(size_t) << " bytes" << std::endl;
            std::cout << "#  * sizeof(size_t*): " << sizeof(size_t*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(numSeqs_t): " << sizeof(numSeqs_t) << " bytes" << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
            std::cout << "#  * strings_: " << (sizeof(S) * capacity_strings_) << " bytes" << std::endl;
            std::cout << "#  * string_cursors_: " << (sizeof(size_t) * capacity_) << " bytes" << std::endl;

        }

    private:
        S* strings_ = nullptr; // overall sequences array, strings do not end with a \0-equivalent
        size_t next_string_ = 0; // position at which the next sequence would be inserted (bit position in strings_)
        size_t* string_cursors_ = nullptr; // array of cursors into strings_ marking the start positions of the individual sequences (bit positions)
        size_t capacity_strings_ = 0; // capacity of strings_ (in number of allocated S)

        numSeqs_t size_ = 0; // current number of amplicons in the collection
        numSeqs_t capacity_ = 0; // capacity of the collection (in terms of the number of amplicons)

        /*
         * Helper methods encapsulating the decoding of a single character and a substring, respectively.
         */
        inline char get_char(const size_t pos) const {
            return rev_acgtu_map[(strings_[pos / (sizeof(S) * 8)] >> ((sizeof(S) * 8) - ((pos % (sizeof(S) * 8)) + 2))) & 3];
        }
        inline std::string get_string(const size_t start, const lenSeqs_t len) const {

            std::string tmp(len, ' ');
            for (size_t i = 0, inner = start % (sizeof(S) * 8), outer = start / (sizeof(S) * 8);
                 i < len;
                 i++, inner = (inner + 2) % (sizeof(S) * 8), outer += (inner == 0)) {

                tmp[i] = rev_acgtu_map[(strings_[outer] >> ((sizeof(S) * 8) - (inner + 2))) & 3];

            }

            return tmp;

        }


        /*
         * Wrapper for access to whole stored strings.
         */
        class TwoBitSeqWrapper : public SequenceWrapper {

        public:
            TwoBitSeqWrapper(const TwoBitStringArray<S>* arr, size_t start, lenSeqs_t len) : arr_(arr), start_(start), len_(len) {
                // nothing else to do
            }

            ~TwoBitSeqWrapper() override = default;

            char operator[](size_t i) const override {
                return arr_->get_char(start_ + 2 * i);
            }

            lenSeqs_t length() const override {
                return len_;
            }

            std::string to_string() const override {
                return arr_->get_string(start_, len_);
            }

            void write(std::ostream& os) const override {
                os << arr_->get_string(start_, len_);
            }

            TwoBitSeqWrapper* clone() const override {
                return new TwoBitSeqWrapper(*this);
            }

            std::unique_ptr<SubstringWrapper> substr(lenSeqs_t start, lenSeqs_t len) const override {
                return std::unique_ptr<SubstringWrapper>(new TwoBitSubWrapper(arr_, start_ + 2 * start, len));
            }

            size_t size_in_bytes() const override {
                return sizeof(TwoBitSeqWrapper);
            }

        private:
            // empty strings are never compared
            bool less(const SequenceWrapper& other) const override {

                auto& o = static_cast<const TwoBitSeqWrapper&>(other);

                bool eq = true;
                lenSeqs_t i = 0;
                for (; eq && i < std::min(len_, other.length()); i++) {
                    eq = arr_->get_char(start_ + 2 * i) == o[i];
                }
                i--;

                return (eq && len_ < other.length()) || (arr_->get_char(start_ + 2 * i) < o[i]);

            }

            bool equals(const SequenceWrapper& other) const override {

                auto& o = static_cast<const TwoBitSeqWrapper&>(other);

                bool eq = len_ == other.length();
                for (auto i = 0; eq && i < len_; i++) {
                    eq = arr_->get_char(start_ + 2 * i) == o[i];
                }

                return eq;

            }


            const TwoBitStringArray<S>* arr_; // container storing the sequence
            size_t start_; // start position of the sequence
            lenSeqs_t len_; // length of the sequence

        };

        /*
         * Wrapper for access to substrings.
         */
        class TwoBitSubWrapper : public SubstringWrapper {

        public:
            TwoBitSubWrapper(const TwoBitStringArray<S>* arr, size_t start, lenSeqs_t len) : arr_(arr), start_(start), len_(len) {
                // nothing else to do
            }

            ~TwoBitSubWrapper() override = default;

            char operator[](size_t i) const override {
                return arr_->get_char(start_ + 2 * i);
            }

            lenSeqs_t length() const override {
                return len_;
            }

            std::string to_string() const override {
                return arr_->get_string(start_, len_);
            }

            TwoBitSubWrapper* clone() const override {
                return new TwoBitSubWrapper(*this);
            }

            void shift() override {
                start_ += 2;
            }

            size_t size_in_bytes() const override {
                return sizeof(TwoBitSubWrapper);
            }

        private:
            // only substrings / segments of equal length are compared
            bool equals(const SubstringWrapper& other) const override {

                auto& o = static_cast<const TwoBitSubWrapper&>(other);

                bool eq = true;
                for (auto i = 0; eq && i < len_; i++) {
                    eq = arr_->get_char(start_ + 2 * i) == o[i];
                }

                return eq;

            }

            const TwoBitStringArray<S>* arr_; // container storing the sequence
            size_t start_; // start position of the sequence
            lenSeqs_t len_; // length of the sequence

        };

    };

    /*
     * Variant of TwoBitStringArray: the cursors pointing into the overall strings array are stored as a DACs sequence.
     */
    template<typename S = size_t>
    class TwoBitDacsStringArray : public StringComponent {
    public:
        TwoBitDacsStringArray() = default;

        TwoBitDacsStringArray(const unsigned long long capacity, const unsigned long long capacity_strings, const ChunkLengths& config) {

            capacity_strings_ = (2 * (capacity_strings + 1 - capacity) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); // + 1 for sentinel, - capacity because final \0 are not stored

            strings_ = new S[capacity_strings_];
            for (auto i = 0; i < capacity_strings_; i++) {
                strings_[i] = 0;
            }
            string_cursors_ = DacsCapacityArray<dacs_t, builder_t, numSeqs_t>(config.seqs_chunk, capacity, config.seqs_level,
                                                                              config.seqs_extend, config.seqs_opt);
            next_string_ = 0;

        }

        ~TwoBitDacsStringArray() {
            delete[] strings_;
        }

        TwoBitDacsStringArray(const TwoBitDacsStringArray& other) { // copy constructor

            capacity_strings_ = other.capacity_strings_;
            strings_ = new S[capacity_strings_];
            memcpy(strings_, other.strings_, sizeof(S) * capacity_strings_);
            string_cursors_ = other.string_cursors_;
            next_string_ = other.next_string_;

        }

        TwoBitDacsStringArray(TwoBitDacsStringArray&& other) noexcept { // move constructor

            capacity_strings_ = other.capacity_strings_;
            next_string_ = other.next_string_;
            strings_ = other.strings_; other.strings_ = nullptr;
            string_cursors_ = std::move(other.string_cursors_);

        }

        TwoBitDacsStringArray& operator=(const TwoBitDacsStringArray& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] strings_;

            // copy new resources
            capacity_strings_ = other.capacity_strings_;
            strings_ = new S[capacity_strings_];
            memcpy(strings_, other.strings_, sizeof(S) * capacity_strings_);
            string_cursors_ = other.string_cursors_;
            next_string_ = other.next_string_;

            return *this;

        }

        TwoBitDacsStringArray& operator=(TwoBitDacsStringArray&& other) noexcept { // move assignment operator

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

            return *this;

        }

        TwoBitDacsStringArray* clone() const override { // deep-copy clone method
            return new TwoBitDacsStringArray(*this);
        }

        inline void push_back(const std::string& s) override {

            extend(s);

            string_cursors_.push_back(next_string_);
            for (size_t i = 0, inner = next_string_ % (sizeof(S) * 8), outer = next_string_ / (sizeof(S) * 8);
                 i < s.length();
                 i++, inner = (inner + 2) % (sizeof(S) * 8), outer += (inner == 0)) {

                strings_[outer] |= static_cast<S>(acgtu_map[s[i]] - 1) << ((sizeof(S) * 8) - 2 - inner);

            }
            next_string_ += 2 * s.length();

        }

        inline std::unique_ptr<SequenceWrapper> get(const numSeqs_t i, const lenSeqs_t len) const override {
            return std::unique_ptr<SequenceWrapper>(new TwoBitSeqWrapper(this, string_cursors_.get(i), len));
        }

        std::string get_str(const numSeqs_t i, const lenSeqs_t len) const override {
            return get_string(string_cursors_.get(i), len);
        }

        inline void extend(const std::string& s) override {

            size_t len_s = 2 * s.length();

            if ((next_string_ + len_s + 2) > capacity_strings_ * (sizeof(S) * 8)) {

                size_t new_capa = (std::max(next_string_ + len_s + 2, 2 * capacity_strings_ * (sizeof(S) * 8)) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8);
                S* tmp1 = new S[new_capa];
                memcpy(tmp1, strings_, capacity_strings_ * sizeof(S));
                delete[] strings_;
                strings_ = tmp1;

                capacity_strings_ = new_capa;

            }

        }

        void finalise() override {
            string_cursors_.finalise();
        }

        void sort(const std::vector<numSeqs_t>& permutation) override {
            string_cursors_.sort(permutation);
        }

        size_t size_in_bytes() const override {
            return sizeof(TwoBitDacsStringArray<S>) // TwoBitDacsStringArray itself
                + sizeof(S) * capacity_strings_ // strings_
                + string_cursors_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>); // string_cursors_ (do not count size of DacsCapacityArray object twice)
        }

        void show_memory() const override {

            std::cout << "#  * TwoBitDacsStringArray " << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * sizeof(TwoBitDacsStringArray<S>): " << sizeof(TwoBitDacsStringArray<S>) << " bytes" << std::endl;
            std::cout << "#  * sizeof(S): " << sizeof(S) << " bytes" << std::endl;
            std::cout << "#  * sizeof(S*): " << sizeof(S*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(size_t): " << sizeof(size_t) << " bytes" << std::endl;
            std::cout << "#  * sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>): " << sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) << " bytes" << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
            std::cout << "#  * strings_: " << (sizeof(S) * capacity_strings_) << " bytes" << std::endl;
            std::cout << "#  * string_cursors_: " << string_cursors_.size_in_bytes() << " bytes" << std::endl;

        }

    private:
        typedef SplitDacs<numSeqs_t, size_t> dacs_t;
        typedef SplitDacsDirectBuilder<numSeqs_t, size_t> builder_t;

        S* strings_ = nullptr; // overall sequences array, strings do not end with a \0-equivalent
        size_t next_string_ = 0; // position at which the next sequence would be inserted (bit position in strings_)
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> string_cursors_; // array of cursors into strings_ marking the start positions of the individual sequences (bit positions)
        size_t capacity_strings_ = 0; // capacity of strings_ (in number of allocated S)


        /*
         * Helper methods encapsulating the decoding of a single character and a substring, respectively.
         */
        inline char get_char(const size_t pos) const {
            return rev_acgtu_map[(strings_[pos / (sizeof(S) * 8)] >> ((sizeof(S) * 8) - ((pos % (sizeof(S) * 8)) + 2))) & 3];
        }
        inline std::string get_string(const size_t start, const lenSeqs_t len) const {

            std::string tmp(len, ' ');
            for (size_t i = 0, inner = start % (sizeof(S) * 8), outer = start / (sizeof(S) * 8);
                 i < len;
                 i++, inner = (inner + 2) % (sizeof(S) * 8), outer += (inner == 0)) {

                tmp[i] = rev_acgtu_map[(strings_[outer] >> ((sizeof(S) * 8) - (inner + 2))) & 3];

            }

            return tmp;

        }


        /*
         * Wrapper for access to whole stored strings.
         */
        class TwoBitSeqWrapper : public SequenceWrapper {

        public:
            TwoBitSeqWrapper(const TwoBitDacsStringArray<S>* arr, size_t start, lenSeqs_t len) : arr_(arr), start_(start), len_(len) {
                // nothing else to do
            }

            ~TwoBitSeqWrapper() override = default;

            char operator[](size_t i) const override {
                return arr_->get_char(start_ + 2 * i);
            }

            lenSeqs_t length() const override {
                return len_;
            }

            std::string to_string() const override {
                return arr_->get_string(start_, len_);
            }

            void write(std::ostream& os) const override {
                os << arr_->get_string(start_, len_);
            }

            TwoBitSeqWrapper* clone() const override {
                return new TwoBitSeqWrapper(*this);
            }

            std::unique_ptr<SubstringWrapper> substr(lenSeqs_t start, lenSeqs_t len) const override {
                return std::unique_ptr<SubstringWrapper>(new TwoBitSubWrapper(arr_, start_ + 2 * start, len));
            }

            size_t size_in_bytes() const override {
                return sizeof(TwoBitSeqWrapper);
            }

        private:
            // empty strings are never compared
            bool less(const SequenceWrapper& other) const override {

                auto& o = static_cast<const TwoBitSeqWrapper&>(other);

                bool eq = true;
                lenSeqs_t i = 0;
                for (; eq && i < std::min(len_, other.length()); i++) {
                    eq = arr_->get_char(start_ + 2 * i) == o[i];
                }
                i--;

                return (eq && len_ < other.length()) || (arr_->get_char(start_ + 2 * i) < o[i]);

            }

            bool equals(const SequenceWrapper& other) const override {

                auto& o = static_cast<const TwoBitSeqWrapper&>(other);

                bool eq = len_ == other.length();
                for (auto i = 0; eq && i < len_; i++) {
                    eq = arr_->get_char(start_ + 2 * i) == o[i];
                }

                return eq;

            }


            const TwoBitDacsStringArray<S>* arr_; // container storing the sequence
            size_t start_; // start position of the sequence
            lenSeqs_t len_; // length of the sequence

        };

        /*
         * Wrapper for access to substrings.
         */
        class TwoBitSubWrapper : public SubstringWrapper {

        public:
            TwoBitSubWrapper(const TwoBitDacsStringArray<S>* arr, size_t start, lenSeqs_t len) : arr_(arr), start_(start), len_(len) {
                // nothing else to do
            }

            ~TwoBitSubWrapper() override = default;

            char operator[](size_t i) const override {
                return arr_->get_char(start_ + 2 * i);
            }

            lenSeqs_t length() const override {
                return len_;
            }

            std::string to_string() const override {
                return arr_->get_string(start_, len_);
            }

            TwoBitSubWrapper* clone() const override {
                return new TwoBitSubWrapper(*this);
            }

            void shift() override {
                start_ += 2;
            }

            size_t size_in_bytes() const override {
                return sizeof(TwoBitSubWrapper);
            }

        private:
            // only substrings / segments of equal length are compared
            bool equals(const SubstringWrapper& other) const override {

                auto& o = static_cast<const TwoBitSubWrapper&>(other);

                bool eq = true;
                for (auto i = 0; eq && i < len_; i++) {
                    eq = arr_->get_char(start_ + 2 * i) == o[i];
                }

                return eq;

            }

            const TwoBitDacsStringArray<S>* arr_; // container storing the sequence
            size_t start_; // start position of the sequence
            lenSeqs_t len_; // length of the sequence

        };

    };

}

#endif //GEFAST_STRINGARRAYS_HPP
