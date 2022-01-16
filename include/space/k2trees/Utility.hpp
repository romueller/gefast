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

#ifndef K2TREES_UTILITY_HPP
#define K2TREES_UTILITY_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace GeFaST {

    /*
     * Refined from https://github.com/romueller/k2trees/blob/master/Utility.hpp.
     */

    /*
     * Position in a two-dimensional matrix.
     * Row and column indices (type I) are 0-based.
     */
    template<typename I = size_t>
    struct Position {

        I row = static_cast<I>(0); // row number
        I col = static_cast<I>(0); // column number

        Position() = default;

        Position(I r, I c) : row(r), col(c) {
            // nothing else to do
        }

        explicit Position(const std::pair<I, I>& pos) : row(pos.first), col(pos.second) {
            // nothing else to do
        }

        bool operator==(const Position& other) const {
            return (row == other.row) && (col == other.col);
        }

        bool operator<(const Position& other) const {
            return (row < other.row) || (row == other.row && col < other.col);
        }

    };


    /*
     * Position in a row of a two-dimensional matrix plus an associated weighted / value of type V.
     * Column indices (type I) are 0-based.
     */
    template<typename V, typename I = size_t>
    struct ValueRowPosition {

        I col = static_cast<I>(0); // column number
        V val = static_cast<V>(0); // value

        ValueRowPosition() = default;

        ValueRowPosition(I c, V v) : col(c), val(v) {
            // nothing else to do
        }

        explicit ValueRowPosition(const std::pair<I, V>& vp) : col(vp.first), val(vp.second) {
            // nothing else to do
        }

        bool operator==(const ValueRowPosition& other) {
            return (col == other.col) && (val == other.val);
        }

    };


    /*
     * Position in a two-dimensional matrix plus an associated weight / value of type V.
     * Row and column indices (type I) are 0-based.
     */
    template<typename V, typename I = size_t>
    struct ValuePosition {

        I row = static_cast<I>(0); // row number
        I col = static_cast<I>(0); // column number
        V val = static_cast<V>(0); // value M[row][col]

        ValuePosition() = default;

        ValuePosition(I r, I c, V v) : row(r), col(c), val(v) {
            // nothing else to do
        }

        ValuePosition(const std::pair<I, I>& pos, V v) : row(pos.first), col(pos.second), val(v) {
            // nothing else to do
        }

        ValuePosition(const Position<I>& pos, V v) : row(pos.row), col(pos.col), val(v) {
            // nothing else to do
        }

        bool operator==(const ValuePosition& other) {
            return (row == other.row) && (col == other.col) && (val == other.val);
        }

        bool operator<(const ValuePosition& other) const {
            return (row < other.row) || (row == other.row && col < other.col);
        }

    };


    /*
     * Collection of partition-related information used in K2Tree-implementations such as PartitionedRectangularK2Tree.
     * Help to map overall indices to relative positions in the partitions.
     */
    template<typename S>
    struct PartitionIndices {

        S partition = static_cast<S>(0); // number / index of the partition
        S row = static_cast<S>(0); // relative row number in the partition
        S col = static_cast<S>(0); // relative column number in the partition

        PartitionIndices() = default;

        PartitionIndices(S p, S r, S c) : partition(p), row(r), col(c) {
            // nothing else to do
        }

    };


    /*
     * Parameters handed over in iterative versions of getting all positions in a row.
     * Used to organise the movement through the tree in a breadth-first manner.
     */
    template<typename S>
    struct SubrowInfo {

        S dq = static_cast<S>(0); // relative column number (on this level)
        S z = static_cast<S>(0); // index in (conceptual concatenation of) T and L

        SubrowInfo() = default;

        SubrowInfo(S d, S index) : dq(d), z(index) {
            // nothing else to do
        }

    };


    /*
     * Parameters handed over in iterative versions of setting all values in a column to null (shallow).
     * Used to organise the movement through the tree in a breadth-first manner.
     */
    template<typename S>
    struct SubcolInfo {

        S dp = static_cast<S>(0); // relative row number (on this level)
        S z = static_cast<S>(0); // index in (conceptual concatenation of) T and L

        SubcolInfo() = default;

        SubcolInfo(S d, S index) : dp(d), z(index) {
            // nothing else to do
        }

    };

    /*
     * Parameters handed over in iterative versions of setting all values in a column to null (deep).
     * Used to organise the movement through the tree in a depth-first manner.
     */
    template<typename S>
    struct SubcolInfoDeep {

        S dp = static_cast<S>(0); // relative row number (on this level)
        S dq = static_cast<S>(0); // relative column number
        S node = static_cast<S>(0); // start index of "node" at which z points
        S z = static_cast<S>(0); // index in (conceptual concatenation of) T and L
        S nr = static_cast<S>(0); // number of rows (on this level)
        S nc = static_cast<S>(0); // number of columns (on this level)
        S k = static_cast<S>(0); // arity (on this level)
        S i = static_cast<S>(0); // index of child considered
        bool zero = true; // indicates whether represented part of matrix is empty

        SubcolInfoDeep() = default;

        SubcolInfoDeep(S p, S q, S node_start, S node_offset, S r, S c, S arity, S child, bool ze)
            : dp(p), dq(q), node(node_start), z(node_start + node_offset), nr(r), nc(c), k(arity), i(child), zero(ze) {
            // nothing else to do
        }

    };


    /*
     * Parameters handed over in iterative versions of getting first positions in a row.
     * Used to organise the movement through the tree in a depth-first manner.
     */
    template<typename S>
    struct ExtendedSubrowInfo {

        S nr = static_cast<S>(0); // number of rows (on this level)
        S nc = static_cast<S>(0); // number of columns (on this level)
        S p = static_cast<S>(0); // relative row number (on this level)
        S dq = static_cast<S>(0); // relative column number (on this level)
        S z = static_cast<S>(0); // index in (conceptual concatenation of) T and L
        S j = static_cast<S>(0); // "child number" (maximum depends on arity)

        ExtendedSubrowInfo() = default;

        ExtendedSubrowInfo(S num_rows, S num_cols, S rel_row, S rel_col, S index, S child)
            : nr(num_rows), nc(num_cols), p(rel_row), dq(rel_col), z(index), j(child) {
            // nothing else to do
        }

    };


    /*
     * Representation of a subproblem / submatrix during the construction of a k^2-tree.
     * Used when constructing the k^2-tree from a list of relation pairs.
     * Describes the part of the matrix and the corresponding range of the input.
     */
    template<typename S = size_t>
    struct Subproblem {

        S first_row = static_cast<S>(0); // first row (x1), inclusive
        S last_row = static_cast<S>(0); // last row (x2), inclusive
        S first_col = static_cast<S>(0); // first column (y1), inclusive
        S last_col = static_cast<S>(0); // last column (y2), inclusive

        S left = static_cast<S>(0); // left border of the elements belonging to the subproblem (a, inclusive index into the input)
        S right = static_cast<S>(0); // right border of the elements belonging to the subproblem (b, exclusive index into the input)

        Subproblem() = default;

        Subproblem(S fr, S lr, S fc, S lc, S l, S r) : first_row(fr), last_row(lr), first_col(fc), last_col(lc), left(l), right(r) {
            // nothing else to do
        }

    };


    /*
     * Helper method for computation of integer log_k(n).
     */
    template<typename S>
    S log_k(const S n, const S k) {

        S pow = 0;
        for (S val = 1; val < n; val *= k, pow++) { }
        return pow;

    }

    /*
     * Helper methods for checking whether all elements of a vector / array have a certain value.
     */
    template<typename T>
    bool is_all(const std::vector<T>& v, const T val) {
        return std::all_of(v.begin(), v.end(), [&val](const T& elem) {return val == elem;});
    }
    template<typename T>
    bool is_all(T* v, const size_t len, const T val) {
        return std::all_of(v, v + len, [&val](const T& elem) {return val == elem;});
    }

    /*
     * Helper methods for checking whether all bits are false.
     */
    bool is_all_zero(const std::vector<bool>& v);
    template<typename B>
    bool is_all_zero(const B& v) { // expects that the container of type B offers operator[]

        bool all_zero = true;
        for (auto i = 0; all_zero && i < v.size(); i++) all_zero = !v[i];
        return all_zero;

    }

    /*
     * Helper method for printing contents of a rank data structure.
     * Expects a rank(index) method with 1-based indices.
     */
    template<typename R>
    void print_ranks(const R& r) {

        for (auto i = 0; i < r.size(); i++) std::cout << r.rank(i + 1) << " " << std::flush;
        std::cout << std::endl;

    }


    /*
     * Dynamic, but naive rank data structure for intermediate steps.
     * Uses 1-based indices. Index 0 is reserved for a guardian value to support rank(0) = 0 without case distinction.
     * Later changes to the bit vector are not seen automatically by the data structure.
     */
    template<typename B, typename S = size_t>
    class NaiveDynamicRank {

    public:
        NaiveDynamicRank() = default;

        explicit NaiveDynamicRank(const B& arr) {

            ranks_ = std::vector<S>(arr.size() + 1);
            ranks_[0] = 0;
            for (auto i = 0; i < arr.size(); i++) {
                ranks_[i + 1] = ranks_[i] + arr[i];
            }

        }

        // return the rank of the 1-based position
        S rank(S pos) const {
            return ranks_[pos];
        }

        // alternative rank method that handles overly large queries
        // as if the bit vector was indefinitely extended by zeros
        S rank_safe(S pos) const {
            return ranks_[std::min(pos, ranks_.size() - 1)];
        }

        /*
         * The next three operations allow to adapt the rank data structure to certain
         * changes of the corresponding bit vector B (0-based indices).
         */

        // operation on bit vector B: change B[pos - 1] from 0 to 1
        void increase_from(S pos) {

            for (auto i = pos; i < ranks_.size(); i++) {
                ranks_[i]++;
            }

        }

        // operation on bit vector B: change B[pos - 1] from 1 to 0
        void decrease_from(S pos) {

            for (auto i = pos; i < ranks_.size(); i++) {
                ranks_[i]--;
            }

        }

        // operation on bit vector B: insert 'num' zeros at position 'pos - 1' in B
        void insert(S pos, S num) {
            ranks_.insert(ranks_.begin() + pos, num, (pos > 0) * ranks_[pos - 1]);
        }

        // length of the underlying bit vector B
        size_t size() const {
            return ranks_.size() - 1;
        }

        // determines memory consumption (in bytes)
        size_t size_in_bytes() const {
            return sizeof(NaiveDynamicRank) // NaiveDynamicRank itself
                + ranks_.capacity() * sizeof(S); // ranks_
        }

    private:
        std::vector<S> ranks_ = {0}; // rank values including the guardian at index 0

    };

    template<typename B, typename S = size_t>
    class RunEncodingRank {

    public:
        RunEncodingRank() = default;

        explicit RunEncodingRank(const B& arr) {

            for (auto i = 0; i < arr.size(); i++) {
                if (arr[i]) {
                    ranks_.push_back(1);
                } else {
                    ranks_.back()++;
                }
            }

        }

        // return the rank of the 1-based position
        S rank(S pos) const {

            S sum = ranks_[0];
            S i = 0;

            while (i < ranks_.size() && sum < pos) {

                i++;
                sum += ranks_[i];

            }

            return i;

        }

        // alternative rank method that handles overly large queries
        // as if the bit vector was indefinitely extended by zeros
        S rank_safe(S pos) const {
            return rank(pos);
        }

        /*
         * The next three operations allow to adapt the rank data structure to certain
         * changes of the corresponding bit vector B (0-based indices).
         */

        // operation on bit vector B: change B[pos - 1] from 0 to 1
        void increase_from(S pos) {

            S sum = ranks_[0];
            S i = 0;
            while (i < ranks_.size() && sum < pos) {

                i++;
                sum += ranks_[i];

            }

            S tmp = sum - pos + 1; // number of values with rank i that are increased

            ranks_[i] -= tmp; // modify count of rank i
            for (auto j = i + 1; j < ranks_.size(); j++) { // integrate increased portion of rank i ...
                std::swap(tmp, ranks_[j]); // ... by pushing back the counts of higher ranks one position
            }
            ranks_.push_back(tmp); // reintegrate the count of the highest rank

        }

        // operation on bit vector B: change B[pos - 1] from 1 to 0
        void decrease_from(S pos) {

            S sum = ranks_[0];
            S i = 0;
            while (i < ranks_.size() && sum < pos) {

                i++;
                sum += ranks_[i];

            }

            ranks_[i - 1] += ranks_[i]; // block of rank i merges with previous block
            for (auto j = i; j < ranks_.size() - 1; j++) {
                ranks_[j] = ranks_[j + 1]; // following block move one position to the front
            }
            ranks_.pop_back(); // previous position of last block is freed

        }

        // operation on bit vector B: insert 'num' zeros at position 'pos - 1' in B
        void insert(S pos, S num) {

            S sum = ranks_[0];
            S i = 0;
            while (i < ranks_.size() && sum < pos - 1) {

                i++;
                sum += ranks_[i];

            }

            ranks_[i] += num; // enlarge the corresponding block

        }

        // length of the underlying bit vector B
        size_t size() const {

            size_t sum = 0;
            for (auto v : ranks_) sum += v;
            return sum;

        }

        // determines memory consumption (in bytes)
        size_t size_in_bytes() const {
            return sizeof(RunEncodingRank) // RunEncodingRank itself
                + ranks_.capacity() * sizeof(S); // ranks_
        }

    private:
        std::vector<S> ranks_ = {0}; // rank values including the guardian at index 0

    };


    /*
     * Dynamic, but naive implementation of a bit vector.
     * Uses 0-based indices and provides some specific dynamic operations.
     * Other that that it is only a wrapper around std::vector<bool>.
     */
    class NaiveBitVector {

    public:
        // create a bit vector of length 'size' with all bits set to false
        explicit NaiveBitVector(size_t size = 0);

        // create a bit vector of length 'size' with all bits set to 'value'
        NaiveBitVector(size_t size, bool value);

        // create a bit vector from std::vector<bool>
        explicit NaiveBitVector(const std::vector<bool>& vec);

        // returns the length of the bit vector in bits
        size_t size() const;

        // checks whether the length is zero
        bool empty() const;

        // sets the i-th bit to 'val'
        void set(size_t i, bool val);

        // returns the value of the i-th bit
        bool get(size_t i) const;

        // returns the value of the i-th bit
        bool operator[](size_t i) const;

        // appends the bits of the provided bit vector at the end
        void append(const std::vector<bool>& bits);

        // overwrites the bits of the bit vector (starting at position i) with the bits
        // of the provided bit vector (does not check whether the bounds)
        size_t set(size_t i, const std::vector<bool>& bits);

        // checks whether all bits with indices i to i + len - 1 are false
        bool is_zero(size_t i, size_t len);

        // determines memory consumption (in bytes)
        size_t size_in_bytes() const;

    private:
        std::vector<bool> bits_; // actual bit vector
    };


}

#endif //K2TREES_UTILITY_HPP
