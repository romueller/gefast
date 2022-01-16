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

#ifndef K2TREES_NAIVEK2TREE_HPP
#define K2TREES_NAIVEK2TREE_HPP

#include "K2Tree.hpp"
#include "Utility.hpp"
#include "../../modes/SpaceLevenshteinMode.hpp"

namespace GeFaST {

    /*
     * Naive implementation of a relation matrix with a K2Tree interface for very small relations.
     *
     * Simply contains a list of the relation pairs.
     *
     * Bool specialisation, which makes use of some simplifications since the only non-null value is 1 / true.
     *
     * Refined class MiniK2Tree<bool> (https://github.com/romueller/k2trees/blob/master/StaticMiniK2Tree.hpp).
     */
    template<typename I>
    class NaiveK2Tree : public virtual K2Tree<bool, I> {

    public:
        typedef bool value_t;
        typedef I index_t;

        typedef typename K2Tree<value_t, index_t>::values_t values_t;
        typedef typename K2Tree<value_t, index_t>::indices_t indices_t;
        typedef typename K2Tree<value_t, index_t>::matrix_t matrix_t;
        typedef typename K2Tree<value_t, index_t>::value_matrix_t value_matrix_t;
        typedef typename K2Tree<value_t, index_t>::row_t row_t;
        typedef typename K2Tree<value_t, index_t>::lists_t lists_t;
        typedef typename K2Tree<value_t, index_t>::value_row_t value_row_t;
        typedef typename K2Tree<value_t, index_t>::value_lists_t value_lists_t;
        typedef typename K2Tree<value_t, index_t>::positions_t positions_t;
        typedef typename K2Tree<value_t, index_t>::value_positions_t value_positions_t;

        typedef typename row_t::const_iterator row_iterator_t;
        typedef typename value_row_t::const_iterator value_row_iterator_t;
        typedef typename positions_t::value_type pos_t;
        typedef typename value_positions_t::value_type value_pos_t;

        NaiveK2Tree() {

            positions_ = nullptr;
            capacity_ = 0;
            num_rows_ = 0;
            num_cols_ = 0;
            size_ = 0;

        }

        /*
         * List-of-pairs-based constructor
         */
        NaiveK2Tree(positions_t& pairs, const size_t num_rows, const size_t num_cols, const size_t kr, const size_t kc) {

            capacity_ = pairs.size();
            size_ = capacity_;
            positions_ = new Position<index_t>[capacity_];

            size_t pos = 0;
            index_t max_row = 0;
            index_t max_col = 0;
            for (auto& p : pairs) {

                positions_[pos] = p;
                pos++;

                max_row = std::max(max_row, p.row);
                max_col = std::max(max_col, p.col);

            }

            size_t h = std::max({static_cast<size_t>(1), log_k(std::max(max_row + 1, num_rows), kr), log_k(std::max(max_col + 1, num_cols), kc)});
            num_rows_ = static_cast<size_t>(std::pow(kr, h));
            num_cols_ = static_cast<size_t>(std::pow(kc, h));

        }

        NaiveK2Tree(positions_t& pairs, const size_t kr, const size_t kc) : NaiveK2Tree(pairs, 0, 0, kr, kc) {
            // nothing else to do
        }

        NaiveK2Tree(positions_t& pairs, const K2TreeParameters& config) : NaiveK2Tree(pairs, 0, 0, config.kr, config.kc) {
            // nothing else to do
        }

        /*
         * List-of-pairs-based constructor similar to the one above, but only a part of the relation matrix is used.
         *
         * Internally, only positions relative to row x and column y are used.
         */
        NaiveK2Tree(positions_t& pairs, const size_t x, const size_t y, const size_t nr, const size_t nc, const size_t l, const size_t r, const size_t kr, const size_t kc) {

            capacity_ = r - l;
            size_ = capacity_;
            positions_ = new Position<index_t>[capacity_];

            for (size_t i = l, pos = 0; i < r; i++, pos++) {

                positions_[pos] = pairs[i];
                positions_[pos].row -= x;
                positions_[pos].col -= y;

            }

            size_t h = std::max({static_cast<size_t>(1), log_k(nr, kr), log_k(nc, kc)});
            num_rows_ = static_cast<size_t>(std::pow(kr, h));
            num_cols_ = static_cast<size_t>(std::pow(kc, h));

        }

        ~NaiveK2Tree() {
            delete[] positions_;
        }

        NaiveK2Tree(const NaiveK2Tree& other) { // copy constructor

            capacity_ = other.capacity_;
            size_ = other.size_;
            num_rows_ = other.num_rows_;
            num_cols_ = other.num_cols_;
            size_ = other.size_;
            positions_ = new Position<index_t>[capacity_];
            for (size_t k = 0; k < capacity_; k++) {
                positions_[k] = other.positions_[k];
            }

        }

        NaiveK2Tree(NaiveK2Tree&& other) noexcept { // move constructor

            capacity_ = other.capacity_;
            size_ = other.size_;
            num_rows_ = other.num_rows_;
            num_cols_ = other.num_cols_;
            positions_ = other.positions_; other.positions_ = nullptr;

        }

        NaiveK2Tree& operator=(const NaiveK2Tree& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] positions_;

            // copy new resources
            capacity_ = other.capacity_;
            size_ = other.size_;
            num_rows_ = other.num_rows_;
            num_cols_ = other.num_cols_;
            positions_ = new Position<index_t>[capacity_];
            for (size_t k = 0; k < capacity_; k++) {
                positions_[k] = other.positions_[k];
            }

            return *this;

        }

        NaiveK2Tree& operator=(NaiveK2Tree&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] positions_;

            // copy / transfer new resources
            capacity_ = other.capacity_;
            size_ = other.size_;
            num_rows_ = other.num_rows_;
            num_cols_ = other.num_cols_;
            positions_ = other.positions_; other.positions_ = nullptr;

            return *this;

        }

        NaiveK2Tree* clone() const override {
            return new NaiveK2Tree(*this);
        }

        size_t get_num_rows() const override {
            return num_rows_;
        }

        size_t get_num_cols() const override {
            return num_cols_;
        }


        bool is_not_null(size_t i, size_t j) override {
            return std::find_if(positions_, positions_ + size_,
                                [i, j](const Position<index_t>& val) {
                                    return val.row == i && val.col == j;
                                }) != positions_ + size_;
        }

        value_t get_value(size_t i, size_t j) override {
            return std::find_if(positions_, positions_ + size_,
                                [i, j](const Position<index_t>& val) {
                                    return val.row == i && val.col == j;
                                }) != positions_ + size_;
        }

        values_t get_successor_values(size_t i) override {

            values_t succs;
            for (size_t k = 0; k < size_; k++) {
                if (positions_[k].row == i) {
                    succs.push_back(true);
                }
            }

            return succs;

        }

        indices_t get_successor_positions(size_t i) override {

            indices_t succs;
            for (size_t k = 0; k < size_; k++) {
                if (positions_[k].row == i) {
                    succs.push_back(positions_[k].col);
                }
            }

            return succs;

        }

        value_positions_t get_successor_value_positions(size_t i) override {

            value_positions_t succs;
            for (size_t k = 0; k < size_; k++) {
                if (positions_[k].row == i) {
                    succs.emplace_back(positions_[k], true);
                }
            }

            return succs;

        }

        values_t get_predecessor_values(size_t j) override {

            values_t preds;
            for (size_t k = 0; k < size_; k++) {
                if (positions_[k].col == j) {
                    preds.push_back(true);
                }
            }

            return preds;

        }

        indices_t get_predecessor_positions(size_t j) override {

            indices_t preds;
            for (size_t k = 0; k < size_; k++) {
                if (positions_[k].col == j) {
                    preds.push_back(positions_[k].row);
                }
            }

            return preds;

        }

        value_positions_t get_predecessor_value_positions(size_t j) override {

            value_positions_t preds;
            for (size_t k = 0; k < size_; k++) {
                if (positions_[k].col == j) {
                    preds.emplace_back(positions_[k], true);
                }
            }

            return preds;

        }

        values_t get_values_in_range(size_t i1, size_t i2, size_t j1, size_t j2) override {

            values_t values;
            for (size_t k = 0; k < size_; k++) {
                if (i1 <= positions_[k].row && positions_[k].row <= i2 && j1 <= positions_[k].col && positions_[k].col <= j2) {
                    values.push_back(true);
                }
            }

            return values;

        }

        positions_t get_positions_in_range(size_t i1, size_t i2, size_t j1, size_t j2) override {

            positions_t pairs;
            for (size_t k = 0; k < size_; k++) {
                if (i1 <= positions_[k].row && positions_[k].row <= i2 && j1 <= positions_[k].col && positions_[k].col <= j2) {
                    pairs.push_back(positions_[k]);
                }
            }

            return pairs;

        }

        value_positions_t get_value_positions_in_range(size_t i1, size_t i2, size_t j1, size_t j2) override {

            value_positions_t pairs;
            for (size_t k = 0; k < size_; k++) {
                if (i1 <= positions_[k].row && positions_[k].row <= i2 && j1 <= positions_[k].col && positions_[k].col <= j2) {
                    pairs.emplace_back(positions_[k], true);
                }
            }

            return pairs;

        }

        values_t get_all_values() override {
            return values_t(size_, true);
        }

        positions_t get_all_positions() override {
            return positions_t(positions_, positions_ + size_);
        }

        value_positions_t get_all_value_positions() override {

            value_positions_t pairs;
            for (size_t k = 0; k < size_; k++) {
                pairs.emplace_back(positions_[k], true);
            }

            return pairs;

        }

        bool contains_value(size_t i1, size_t i2, size_t j1, size_t j2) override {

            bool flag = false;
            for (size_t k = 0; k < size_ && !flag; k++) {
                flag = i1 <= positions_[k].row && positions_[k].row <= i2 && j1 <= positions_[k].col && positions_[k].col <= j2;
            }

            return flag;

        }

        size_t count_values() override {
            return size_;
        }


        size_t size_in_bytes() const override {
            return sizeof(NaiveK2Tree<I>) // NaiveK2Tree itself
                + capacity_ * sizeof(Position<index_t>); // positions_
        }

        void print(bool all) override {

            std::cout << "### Parameters ###" << std::endl;
            std::cout << "num_rows  = " << get_num_rows() << std::endl;
            std::cout << "num_cols  = " << get_num_cols() << std::endl;

            if (all) {

                std::cout << "### Positions & values ###" << std::endl;
                for (size_t k = 0; k < size_; k++) {
                    std::cout << "(" << positions_[k].row << ", " << positions_[k].col << ", " << true << ") ";
                }
                std::cout << std::endl << std::endl;

            }

        }

        void set_null(size_t i, size_t j, bool deep = true) override {

            auto iter = std::find_if(positions_, positions_ + size_,
                                     [i, j](const Position<index_t>& val) {
                                         return val.row == i && val.col == j;
                                     });

            if (iter != positions_ + size_) {

                positions_[iter - positions_] = positions_[size_ - 1];
                size_--;

            }

        }

        void set_null_deep(size_t i, size_t j) override {
            set_null(i, j);
        }

        void set_null_shallow(size_t i, size_t j) override {
            set_null(i, j);
        }

        void set_null_column(size_t j, bool = true) override {


            auto iter = std::find_if(positions_, positions_ + size_,
                                     [j](const Position<index_t>& val) {
                                         return val.col == j;
                                     });

            while (iter != positions_ + size_) {

                positions_[iter - positions_] = positions_[size_ - 1];
                size_--;

                iter = std::find_if(iter, positions_ + size_,
                                    [j](const Position<index_t>& val) {
                                        return val.col == j;
                                    });

            }

        }

        void set_null_column_deep(size_t j) override {
            set_null_column(j);
        }

        void set_null_column_shallow(size_t j) override {
            set_null_column(j);
        }

        index_t get_first_successor(index_t i) override {

            index_t min = get_num_cols();
            for (size_t k = 0; k < size_; k++) {
                if (positions_[k].row == i) {
                    min = std::min(min, positions_[k].col);
                }
            }

            return min;

        }


        /*
         * Method aliases using "relation nomenclature" (similar to the names proposed by Brisaboa et al.)
         */

        bool are_related(index_t i, index_t j) override {
            return is_not_null(i, j);
        }

        indices_t get_successors(index_t i) override {
            return get_successor_positions(i);
        }

        indices_t get_predecessors(index_t j) override {
            return get_predecessor_positions(j);
        }

        positions_t get_range(index_t i1, index_t i2, index_t j1, index_t j2) override {
            return get_positions_in_range(i1, i2, j1, j2);
        }

        bool contains_link(index_t i1, index_t i2, index_t j1, index_t j2) override {
            return contains_value(i1, i2, j1, j2);
        }

        size_t count_links() override {
            return count_values();
        }



    private:
        Position<index_t>* positions_; // positions of all relation pairs
        size_t capacity_; // initial number of relation pairs

        size_t num_rows_; // number of rows in the represented relation matrix
        size_t num_cols_; // number of columns in the represented relation matrix
        size_t size_; // number of relation pairs

    };

}

#endif //K2TREES_NAIVEK2TREE_HPP
