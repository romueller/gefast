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

#ifndef K2TREES_PARTITIONEDRECTANGULARK2TREE_HPP
#define K2TREES_PARTITIONEDRECTANGULARK2TREE_HPP

#include <queue>

#include "K2Tree.hpp"
#include "RectangularK2Tree.hpp"
#include "Utility.hpp"
#include "../../modes/SpaceLevenshteinMode.hpp"

namespace GeFaST {

    /*
     * Uneven rectangular implementation of K2Tree.
     *
     * Uses two different arities for rows and columns (kr, kc) and allows for a different
     * number of levels over rows and columns. This effectively leads to a partitioning
     * of the whole relation matrix into several RectangularK2Tree instances.
     * The described relation matrix is rectangular with edge lengths of num_rows and num_cols,
     * where num_rows (num_cols) is the smallest power of kr (kc) that exceeds the row (column) numbers
     * of all relation pairs.
     *
     * Bool specialisation, which makes use of some simplifications since the only non-null value is 1 / true.
     *
     * Refined class UnevenKrKcTree<bool> (https://github.com/romueller/k2trees/blob/master/StaticUnevenRectangularTree.hpp).
     */
    template<typename B, typename R, typename I>
    class PartitionedRectangularK2Tree : public virtual K2Tree<bool, I> {

    public:
        typedef bool value_t;
        typedef B bit_vector_t;
        typedef R rank_t;
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

        PartitionedRectangularK2Tree() {

            kr_ = 0;
            kc_ = 0;
            hr_ = 0;
            hc_ = 0;
            num_rows_ = 0;
            num_cols_ = 0;
            partitions_ = nullptr;
            partition_size_ = 0;
            num_partitions_ = 0;

        }

        /*
         * List-of-pairs-based constructor (based on section 3.3.5. of Brisaboa et al.)
         */
        PartitionedRectangularK2Tree(positions_t& pairs, const size_t num_rows, const size_t num_cols, const size_t kr, const size_t kc, const bool perm_rank) {

            size_t max_row = 0;
            size_t max_col = 0;
            for (auto& p : pairs) {
                max_row = std::max(max_row, p.row);
                max_col = std::max(max_col, p.col);
            }

            kr_ = kr;
            kc_ = kc;
            hr_ = std::max(static_cast<size_t>(1), log_k(std::max(max_row + 1, num_rows), kr_));
            hc_ = std::max(static_cast<size_t>(1), log_k(std::max(max_col + 1, num_cols), kc_));
            num_rows_ = static_cast<size_t>(std::pow(kr_, hr_));
            num_cols_ = static_cast<size_t>(std::pow(kc_, hc_));

            if (hc_ > hr_) {

                partition_size_ = static_cast<size_t>(std::pow(kc_, hr_));
                num_partitions_ = num_cols_ / partition_size_;
                partitions_ = new RectangularK2Tree<bit_vector_t, rank_t, index_t>*[num_partitions_];

                Subproblem<size_t> sp(0, num_rows_ - 1, 0, num_cols_ - 1, 0, pairs.size());
                std::vector<std::pair<size_t, size_t>> intervals(num_partitions_);
                counting_sort(pairs, intervals, sp, num_rows_, partition_size_, num_partitions_);

                for (size_t i = 0; i < num_partitions_; i++) {
                    partitions_[i] = new RectangularK2Tree<bit_vector_t, rank_t, index_t>(pairs, 0, i * partition_size_, num_rows_, partition_size_, intervals[i].first, intervals[i].second, kr_, kc_, perm_rank);
                }


            } else {

                partition_size_ = static_cast<size_t>(std::pow(kr_, hc_));
                num_partitions_ = num_rows_ / partition_size_;
                partitions_ = new RectangularK2Tree<bit_vector_t, rank_t, index_t>*[num_partitions_];

                Subproblem<size_t> sp(0, num_rows_ - 1, 0, num_cols_ - 1, 0, pairs.size());
                std::vector<std::pair<size_t, size_t>> intervals(num_partitions_);
                counting_sort(pairs, intervals, sp, partition_size_, num_cols_, num_partitions_);

                for (size_t j = 0; j < num_partitions_; j++) {
                    partitions_[j] = new RectangularK2Tree<bit_vector_t, rank_t, index_t>(pairs, j * partition_size_, 0, partition_size_, num_cols_, intervals[j].first, intervals[j].second, kr_, kc_, perm_rank);
                }

            }

            for (size_t k = 0; k < num_partitions_; k++) {

                if (partitions_[k]->get_num_rows() == 0) {

                    delete partitions_[k];
                    partitions_[k] = nullptr;

                }

            }

        }

        PartitionedRectangularK2Tree(positions_t& pairs, const size_t kr, const size_t kc, const bool perm_rank) : PartitionedRectangularK2Tree(pairs, 0, 0, kr, kc, perm_rank) {
            // nothing else to do
        }

        PartitionedRectangularK2Tree(positions_t& pairs, const K2TreeParameters& config) : PartitionedRectangularK2Tree(pairs, 0, 0, config.kr, config.kc, config.deep) {
            // nothing else to do
        }

        ~PartitionedRectangularK2Tree() {

            for (size_t k = 0; k < num_partitions_; k++) {
                delete partitions_[k];
            }
            delete[] partitions_;

        }

        PartitionedRectangularK2Tree(const PartitionedRectangularK2Tree& other) { // copy constructor

            hr_ = other.hr_;
            hc_ = other.hc_;
            kr_ = other.kr_;
            kc_ = other.kc_;
            num_rows_ = other.num_rows_;
            num_cols_ = other.num_cols_;
            null_ = other.null_;
            partition_size_ = other.partition_size_;
            num_partitions_ = other.num_partitions_;

            partitions_ = new RectangularK2Tree<bit_vector_t, rank_t, index_t>*[num_partitions_];
            for (size_t k = 0; k < num_partitions_; k++) {
                partitions_[k] = new RectangularK2Tree<bit_vector_t, rank_t, index_t>(*other.partitions_[k]);
            }

        }

        PartitionedRectangularK2Tree(PartitionedRectangularK2Tree&& other) noexcept { // move constructor

            hr_ = other.hr_;
            hc_ = other.hc_;
            kr_ = other.kr_;
            kc_ = other.kc_;
            num_rows_ = other.num_rows_;
            num_cols_ = other.num_cols_;
            null_ = other.null_;
            partition_size_ = other.partition_size_;
            num_partitions_ = other.num_partitions_;

            partitions_ = other.partitions_; other.partitions_ = nullptr;
            other.num_partitions_ = 0;

        }

        PartitionedRectangularK2Tree& operator=(const PartitionedRectangularK2Tree& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            for (size_t k = 0; k < num_partitions_; k++) {
                delete partitions_[k];
            }
            delete[] partitions_;

            // copy new resources
            hr_ = other.hr_;
            hc_ = other.hc_;
            kr_ = other.kr_;
            kc_ = other.kc_;
            num_rows_ = other.num_rows_;
            num_cols_ = other.num_cols_;
            null_ = other.null_;
            partition_size_ = other.partition_size_;
            num_partitions_ = other.num_partitions_;

            partitions_ = new RectangularK2Tree<bit_vector_t, rank_t, index_t>*[num_partitions_];
            for (size_t k = 0; k < num_partitions_; k++) {
                partitions_[k] = new RectangularK2Tree<bit_vector_t, rank_t, index_t>(*other.partitions_[k]);
            }

            return *this;

        }

        PartitionedRectangularK2Tree& operator=(PartitionedRectangularK2Tree&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            for (size_t k = 0; k < num_partitions_; k++) {
                delete partitions_[k];
            }
            delete[] partitions_;

            // copy / transfer new resources
            hr_ = other.hr_;
            hc_ = other.hc_;
            kr_ = other.kr_;
            kc_ = other.kc_;
            num_rows_ = other.num_rows_;
            num_cols_ = other.num_cols_;
            null_ = other.null_;
            partition_size_ = other.partition_size_;
            num_partitions_ = other.num_partitions_;

            partitions_ = other.partitions_; other.partitions_ = nullptr;
            other.num_partitions_ = 0;

            return *this;

        }

        PartitionedRectangularK2Tree* clone() const override {
            return new PartitionedRectangularK2Tree(*this);
        }


        // returns the row height of the K2Tree
        size_t get_row_height() const {
            return hr_;
        }

        // returns the column height of the K2Tree
        size_t get_col_height() const {
            return hc_;
        }

        // returns the row arity of the K2Tree
        size_t get_row_arity() const {
            return kr_;
        }

        // returns the column arity of the K2Tree
        size_t get_col_arity() const {
            return kc_;
        }

        size_t get_num_rows() const override {
            return num_rows_;
        }

        size_t get_num_cols() const override {
            return num_cols_;
        }


        bool are_related(index_t i, index_t j) override {

            auto pis = determine_indices(i, j);
            auto p = partitions_[pis.partition];

            return (p != nullptr) && p->are_related(pis.row, pis.col);

        }

        indices_t get_successors(index_t i) override {

            indices_t succs;

            if (hc_ > hr_) {

                size_t offset = 0;
                for (size_t k = 0; k < num_partitions_; k++, offset += partition_size_) {

                    auto p = partitions_[k];
                    if (p != nullptr) {

                        auto tmp = p->get_successors(i);
                        for (size_t l = 0; l < tmp.size(); l++) {
                            tmp[l] += offset;
                        }

                        succs.reserve(succs.size() + tmp.size());
                        std::move(tmp.begin(), tmp.end(), std::back_inserter(succs));

                    }

                }

            } else {

                auto pis = determine_indices(i, 0);
                auto p = partitions_[pis.partition];
                if (p != nullptr) {
                    succs = p->get_successors(pis.row);
                }

            }

            return succs;

        }

        indices_t get_predecessors(index_t j) override {

            indices_t preds;

            if (hc_ < hr_) {

                size_t offset = 0;
                for (size_t k = 0; k < num_partitions_; k++, offset += partition_size_) {

                    auto p = partitions_[k];
                    if (p != nullptr) {

                        auto tmp = p->get_predecessors(j);
                        for (size_t l = 0; l < tmp.size(); l++) {
                            tmp[l] += offset;
                        }

                        preds.reserve(preds.size() + tmp.size());
                        std::move(tmp.begin(), tmp.end(), std::back_inserter(preds));

                    }

                }

            } else {

                auto pis = determine_indices(0, j);
                auto p = partitions_[pis.partition];
                if (p != nullptr) {
                    preds = p->get_predecessors(pis.col);
                }

            }

            return preds;

        }

        positions_t get_range(index_t i1, index_t i2, index_t j1, index_t j2) override {

            positions_t pairs;

            auto upper_left = determine_indices(i1, j1);
            auto lower_right = determine_indices(i2, j2);

            // range falls completely within one partition
            if (upper_left.partition == lower_right.partition) {

                auto p = partitions_[upper_left.partition];
                if (p != nullptr) {

                    pairs = p->get_range(upper_left.row, lower_right.row, upper_left.col, lower_right.col);

                    size_t offset = partition_size_ * upper_left.partition;

                    if (hc_ > hr_) {

                        for (auto& e : pairs) {
                            e.col += offset;
                        }

                    } else {

                        for (auto& e : pairs) {
                            e.row += offset;
                        }

                    }

                }

                return pairs;

            }

            // range spans multiple partitions
            std::vector<positions_t> tmp;
            tmp.reserve(lower_right.partition - upper_left.partition + 1);
            size_t offset = partition_size_ * upper_left.partition;

            if (hc_ > hr_) {

                // first partition (partially spanned)
                auto p = partitions_[upper_left.partition];
                if (p != nullptr) {

                    tmp.push_back(p->get_range(upper_left.row, lower_right.row, upper_left.col, partition_size_ - 1));
                    for (auto& e : tmp.back()) {
                        e.col += offset;
                    }

                }

                // intermediate partition (fully spanned, if any)
                offset += partition_size_;
                for (size_t k = upper_left.partition + 1; k < lower_right.partition; k++, offset += partition_size_) {

                    p = partitions_[k];
                    if (p != nullptr) {

                        tmp.push_back(p->get_range(upper_left.row, lower_right.row, 0, partition_size_ - 1));
                        for (auto& e : tmp.back()) {
                            e.col += offset;
                        }

                    }

                }

                // last partition (partially spanned)
                p = partitions_[lower_right.partition];
                if (p != nullptr) {

                    tmp.push_back(p->get_range(upper_left.row, lower_right.row, 0, lower_right.col));
                    for (auto& e : tmp.back()) {
                        e.col += offset;
                    }

                }

            } else {

                // first partition (partially spanned)
                auto p = partitions_[upper_left.partition];
                if (p != nullptr) {

                    tmp.push_back(p->get_range(upper_left.row, partition_size_ - 1, upper_left.col, lower_right.col));
                    for (auto& e : tmp.back()) {
                        e.row += offset;
                    }

                }

                // intermediate partition (fully spanned, if any)
                offset += partition_size_;
                for (size_t k = upper_left.partition + 1; k < lower_right.partition; k++, offset += partition_size_) {

                    p = partitions_[k];
                    if (p != nullptr) {

                        tmp.push_back(p->get_range(0, partition_size_ - 1, upper_left.col, lower_right.col));
                        for (auto& e : tmp.back()) {
                            e.row += offset;
                        }

                    }

                }

                // last partition (partially spanned)
                p = partitions_[lower_right.partition];
                if (p != nullptr) {

                    tmp.push_back(p->get_range(0, lower_right.row, upper_left.col, lower_right.col));
                    for (auto& e : tmp.back()) {
                        e.row += offset;
                    }

                }

            }

            // flatten results
            size_t to_reserve = pairs.size();
            for (auto& v : tmp) {
                to_reserve += v.size();
            }
            pairs.reserve(to_reserve);

            for (auto& v : tmp) {

                std::move(v.begin(), v.end(), std::back_inserter(pairs));
                v.clear();
                v.shrink_to_fit();

            }

            return pairs;

        }

        bool contains_link(index_t i1, index_t i2, index_t j1, index_t j2) override {

            auto upper_left = determine_indices(i1, j1);
            auto lower_right = determine_indices(i2, j2);

            // range falls completely within one partition
            if (upper_left.partition == lower_right.partition) {

                auto p = partitions_[upper_left.partition];
                return (p != nullptr) && p->contains_link(upper_left.row, lower_right.row, upper_left.col, lower_right.col);

            }

            // range spans multiple partitions
            bool found;

            if (hc_ > hr_) {

                // first partition (partially spanned)
                auto p = partitions_[upper_left.partition];
                found = (p != nullptr) && p->contains_link(upper_left.row, lower_right.row, upper_left.col, partition_size_ - 1);

                // intermediate partition (fully spanned, if any)
                for (size_t k = upper_left.partition + 1; (k < lower_right.partition) && !found; k++) {

                    p = partitions_[k];
                    found = (p != nullptr) && p->contains_link(upper_left.row, lower_right.row, 0, partition_size_ - 1);

                }

                // last partition (partially spanned)
                p = partitions_[lower_right.partition];
                found = found || ((p != nullptr) && p->contains_link(upper_left.row, lower_right.row, 0, lower_right.col));

            } else {

                // first partition (partially spanned)
                auto p = partitions_[upper_left.partition];
                found =  (p != nullptr) && p->contains_link(upper_left.row, partition_size_ - 1, upper_left.col, lower_right.col);

                // intermediate partition (fully spanned, if any)
                for (size_t k = upper_left.partition + 1; (k < lower_right.partition) && !found; k++) {

                    p = partitions_[k];
                    found = (p != nullptr) && p->contains_link(0, partition_size_ - 1, upper_left.col, lower_right.col);

                }

                // last partition (partially spanned)
                p = partitions_[lower_right.partition];
                found = found ||  ((p != nullptr) && p->contains_link(0, lower_right.row, upper_left.col, lower_right.col));

            }

            return found;

        }

        size_t count_links() override {

            size_t cnt = 0;
            for (size_t k = 0; k < num_partitions_; k++) {

                auto p = partitions_[k];
                cnt += (p != nullptr) ? p->count_links() : 0;

            }

            return cnt;

        }


        /*
         * General methods for completeness' sake (are redundant / useless for bool)
         */

        bool is_not_null(index_t i, index_t j) override {
            return are_related(i, j);
        }

        value_t get_value(index_t i, index_t j) override {
            return are_related(i, j);
        }

        values_t get_successor_values(index_t i) override {
            return values_t(get_successors(i).size(), true);
        }

        indices_t get_successor_positions(index_t i) override {
            return get_successors(i);
        }

        value_positions_t get_successor_value_positions(index_t i) override {

            auto pos = get_successors(i);

            value_positions_t succs;
            for (auto j : pos) {
                succs.emplace_back(i, j, true);
            }

            return succs;

        }

        values_t get_predecessor_values(index_t j) override {
            return values_t(get_predecessors(j).size(), true);
        }

        indices_t get_predecessor_positions(index_t j) override {
            return get_predecessors(j);
        }

        value_positions_t get_predecessor_value_positions(index_t j) override {

            auto pos = get_predecessors(j);

            value_positions_t preds;
            for (auto i : pos) {
                preds.emplace_back(i, j, true);
            }

            return preds;

        }

        values_t get_values_in_range(index_t i1, index_t i2, index_t j1, index_t j2) override {
            return values_t(get_range(i1, i2, j1, j2).size(), true);
        }

        positions_t get_positions_in_range(index_t i1, index_t i2, index_t j1, index_t j2) override {
            return get_range(i1, i2, j1, j2);
        }

        value_positions_t get_value_positions_in_range(index_t i1, index_t i2, index_t j1, index_t j2) override {

            auto pos = get_range(i1, i2, j1, j2);

            value_positions_t pairs;
            for (auto& p : pos) {
                pairs.emplace_back(p.row, p.col, true);
            }

            return pairs;

        }

        values_t get_all_values() override {
            return values_t(count_links(), true);
        }

        positions_t get_all_positions() override {
            return get_range(0, num_rows_ - 1, 0, num_cols_ - 1);
        }

        value_positions_t get_all_value_positions() override {

            auto pos = get_all_positions();

            value_positions_t pairs;
            for (auto& p : pos) {
                pairs.emplace_back(p.row, p.col, true);
            }

            return pairs;

        }

        bool contains_value(index_t i1, index_t i2, index_t j1, index_t j2) override {
            return contains_link(i1, i2, j1, j2);
        }

        size_t count_values() override {
            return count_links();
        }


        size_t size_in_bytes() const override {

            size_t res = sizeof(PartitionedRectangularK2Tree<B, R, I>) // PartitionedRectangularK2Tree itself
                    + num_partitions_ * sizeof(RectangularK2Tree<B, R, I>*); // partitions_ (pointers)
            for (size_t i = 0; i < num_partitions_; i++) {
                res += partitions_[i]->size_in_bytes(); // partitions_ (k^2-trees)
            }

            return res;

        }

        void print(bool all) override {

            std::cout << "### Parameters ###" << std::endl;
            std::cout << "hr  = " << hr_ << std::endl;
            std::cout << "hc  = " << hc_ << std::endl;
            std::cout << "kr  = " << kr_ << std::endl;
            std::cout << "kc  = " << kc_ << std::endl;
            std::cout << "num_rows = " << num_rows_ << std::endl;
            std::cout << "num_cols = " << num_cols_ << std::endl;
            std::cout << "partition_size = " << partition_size_ << std::endl;
            std::cout << "num_partitions = " << num_partitions_ << std::endl;

            if (all) {

                for (size_t k = 0; k < num_partitions_; k++) {

                    std::cout << "===== Partition " << k << " =====" << std::endl;
                    auto p = partitions_[k];
                    if (p != nullptr) {
                        p->print(true);
                    } else {
                        std::cout << "((ALL NULL))" << std::endl;
                    }

                }

            }

        }

        void set_null(size_t i, size_t j, bool deep) override {

            auto pis = determine_indices(i, j);
            auto p = partitions_[pis.partition];

            if (p != nullptr) {
                p->set_null(pis.row, pis.col, deep);
            }

        }

        // note: can "invalidate" the data structure (contains_link() probably won't work correctly afterwards)
        void set_null_deep(size_t i, size_t j) override {

            auto pis = determine_indices(i, j);
            auto p = partitions_[pis.partition];

            if (p != nullptr) {
                p->set_null_deep(pis.row, pis.col);
            }

        }

        // note: can "invalidate" the data structure (contains_link() probably won't work correctly afterwards)
        void set_null_shallow(size_t i, size_t j) override {

            auto pis = determine_indices(i, j);
            auto p = partitions_[pis.partition];

            if (p != nullptr) {
                p->set_null_shallow(pis.row, pis.col);
            }

        }

        void set_null_column(index_t j, bool deep) override {

            if (hc_ < hr_) {

                for (size_t k = 0; k < num_partitions_; k++) {

                    auto p = partitions_[k];
                    if (p != nullptr) {
                        p->set_null_column(j, deep);
                    }

                }

            } else {

                auto pis = determine_indices(0, j);
                auto p = partitions_[pis.partition];
                if (p != nullptr) {
                    p->set_null_column(pis.col, deep);
                }

            }

        }

        // note: can "invalidate" the data structure (contains_link() probably won't work correctly afterwards)
        void set_null_column_deep(index_t j) override {

            if (hc_ < hr_) {

                for (size_t k = 0; k < num_partitions_; k++) {

                    auto p = partitions_[k];
                    if (p != nullptr) {
                        p->set_null_column_deep(j);
                    }

                }

            } else {

                auto pis = determine_indices(0, j);
                auto p = partitions_[pis.partition];
                if (p != nullptr) {
                    p->set_null_column_deep(pis.col);
                }

            }

        }

        // note: can "invalidate" the data structure (contains_link() probably won't work correctly afterwards)
        void set_null_column_shallow(index_t j) override {

            if (hc_ < hr_) {

                for (size_t k = 0; k < num_partitions_; k++) {

                    auto p = partitions_[k];
                    if (p != nullptr) {
                        p->set_null_column_shallow(j);
                    }

                }

            } else {

                auto pis = determine_indices(0, j);
                auto p = partitions_[pis.partition];
                if (p != nullptr) {
                    p->set_null_column_shallow(pis.col);
                }

            }

        }

        size_t get_first_successor(size_t i) override {

            size_t pos = num_cols_;

            if (hc_ > hr_) {

                size_t offset = 0;
                for (size_t k = 0; k < num_partitions_ && pos == num_cols_; k++, offset += partition_size_) {

                    auto p = partitions_[k];
                    if (p != nullptr) {

                        auto tmp = p->get_first_successor(i);
                        if (tmp != p->get_num_cols()) {
                            pos = offset + tmp;
                        }

                    }

                }

            } else {

                auto pis = determine_indices(i, 0);
                auto p = partitions_[pis.partition];
                if (p != nullptr) {
                    pos = p->get_first_successor(pis.row);
                }

            }

            return pos;

        }



    private:
        size_t hr_; // row height of the K2Tree
        size_t hc_; // column height of the K2Tree
        size_t kr_; // row arity of the K2Tree
        size_t kc_; // column arity of the K2Tree
        size_t num_rows_; // number of rows in the represented relation matrix
        size_t num_cols_; // number of columns in the represented relation matrix

        RectangularK2Tree<bit_vector_t, rank_t, index_t>** partitions_; // representations of the partitions / submatrices
        size_t partition_size_; // number of rows (columns) per partition in a vertical (horizontal) partitioning
        size_t num_partitions_; // number of partitions

        value_t null_ = false; // null value


        /* helper methods for mapping (overall) indices to positions in the partitions */

        size_t determine_partition(size_t i, size_t j) {
            return (hc_ > hr_) ? (j / partition_size_) : (i / partition_size_);
        }

        PartitionIndices<size_t> determine_indices(size_t i, size_t j) {
            return (hc_ > hr_) ? PartitionIndices<size_t>(j / partition_size_, i, j % partition_size_) : PartitionIndices<size_t>(i / partition_size_, i % partition_size_, j);
        }

        /* helper methods for inplace construction from single list of pairs */

        size_t compute_key(const pos_t& pair, const Subproblem<size_t>& sp, size_t width_row, size_t width_col) {
            return ((pair.row - sp.first_row) / width_row) + (pair.col - sp.first_col) / width_col;
        }

        void counting_sort(positions_t& pairs, std::vector<std::pair<size_t, size_t>>& intervals, const Subproblem<size_t>& sp, size_t width_row, size_t width_col, size_t sup) {

            std::vector<size_t> counts(sup);

            // determine key frequencies
            for (size_t i = sp.left; i < sp.right; i++) {
                counts[compute_key(pairs[i], sp, width_row, width_col)]++;
            }

            // determine starting index for each key
            size_t total = 0;
            size_t tmp;

            for (size_t key = 0; key < sup; key++) {

                tmp = counts[key];
                counts[key] = total;
                total += tmp;

                intervals[key].first = counts[key];
                intervals[key].second = total;

            }

            // reorder pairs of current subproblem
            positions_t tmp_pairs(sp.right - sp.left + 1);
            for (size_t i = sp.left; i < sp.right; i++) {

                tmp_pairs[counts[compute_key(pairs[i], sp, width_row, width_col)]] = pairs[i];
                counts[compute_key(pairs[i], sp, width_row, width_col)]++;

            }

            for (size_t i = sp.left; i < sp.right; i++) {
                pairs[i] = tmp_pairs[i - sp.left];
            }

        }

    };

}

#endif //K2TREES_PARTITIONEDRECTANGULARK2TREE_HPP
