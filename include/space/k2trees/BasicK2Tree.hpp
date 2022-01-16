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

#ifndef K2TREES_BASICK2TREE_HPP
#define K2TREES_BASICK2TREE_HPP

#include <queue>
#include <stack>

#include "K2Tree.hpp"
#include "Utility.hpp"
#include "../../modes/SpaceLevenshteinMode.hpp"

namespace GeFaST {

    /*
     * Simplest implementation of K2Tree.
     *
     * Uses the same arity (k) on all levels and for both rows and columns.
     * The described relation matrix is quadratic with an edge length of n_prime,
     * where n_prime is the smallest power of k that exceeds the row / column numbers
     * of all relation pairs.
     *
     * Bool specialisation, which makes use of some simplifications since the only non-null value is 1 / true.
     *
     * Refined class BasicK2Tree<bool> (https://github.com/romueller/k2trees/blob/master/StaticBasicTree.hpp).
     */
    template<typename B, typename R, typename I>
    class BasicK2Tree : public K2Tree<bool, I> {

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

        BasicK2Tree() {

            h_ = 0;
            k_ = 0;
            n_prime_ = 0;

        }

        /*
         * List-of-pairs-based constructor (based on section 3.3.5. of Brisaboa et al.)
         */
        BasicK2Tree(positions_t& pairs, const size_t num_rows, const size_t num_cols, const size_t k, const bool perm_rank) {

            index_t max_index = 0;
            for (auto& p : pairs) {
                max_index = std::max({max_index, p.row, p.col});
            }

            k_ = k;
            h_ = std::max(static_cast<size_t>(1), log_k(std::max({max_index + 1, num_rows, num_cols}), k_));
            n_prime_ = static_cast<size_t>(std::pow(k_, h_));

            if (!pairs.empty()) {
                build_from_pairs_inplace(pairs);
            }

            perm_rank_ = perm_rank;
            ranks_ = rank_t(tree_interior_, perm_rank_);

        }

        BasicK2Tree(positions_t& pairs, const size_t k, const bool perm_rank) : BasicK2Tree(pairs, 0, 0, k, perm_rank) {
            // nothing else to do
        }

        BasicK2Tree(positions_t& pairs, const K2TreeParameters& config) : BasicK2Tree(pairs, 0, 0, config.k, config.deep) {
            // nothing else to do
        }

        BasicK2Tree(const BasicK2Tree& other) { // copy constructor

            h_ = other.h_;
            k_ = other.k_;
            n_prime_ = other.n_prime_;
            null_ = other.null_;
            perm_rank_ = other.perm_rank_;

            tree_interior_ = other.tree_interior_;
            leaves_ = other.leaves_;
            ranks_ = rank_t(tree_interior_, perm_rank_);

        }

        BasicK2Tree(BasicK2Tree&& other) noexcept : tree_interior_(std::move(other.tree_interior_)), leaves_(std::move(other.leaves_)), ranks_(std::move(other.ranks_)) { // move constructor

            h_ = other.h_;
            k_ = other.k_;
            n_prime_ = other.n_prime_;
            null_ = other.null_;
            perm_rank_ = other.perm_rank_;

        }

        BasicK2Tree& operator=(const BasicK2Tree& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            h_ = other.h_;
            k_ = other.k_;
            n_prime_ = other.n_prime_;
            null_ = other.null_;
            perm_rank_ = other.perm_rank_;

            tree_interior_ = other.tree_interior_;
            leaves_ = other.leaves_;
            ranks_ = rank_t(tree_interior_, perm_rank_);

            return *this;

        }

        BasicK2Tree& operator=(BasicK2Tree&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            h_ = other.h_;
            k_ = other.k_;
            n_prime_ = other.n_prime_;
            null_ = other.null_;
            perm_rank_ = other.perm_rank_;

            tree_interior_ = std::move(other.tree_interior_);
            leaves_ = std::move(other.leaves_);
            ranks_ = std::move(other.ranks_);

            return *this;

        }

        BasicK2Tree* clone() const override {
            return new BasicK2Tree(*this);
        }


        // returns the height of the K2Tree
        size_t get_height() const {
            return h_;
        }

        // returns the arity of the K2Tree
        size_t get_arity() const {
            return k_;
        }

        index_t get_num_rows() const override {
            return n_prime_;
        }

        index_t get_num_cols() const override {
            return n_prime_;
        }


        bool are_related(index_t i, index_t j) override {
            return check_link_init(i, j);
        }

        indices_t get_successors(index_t i) override {

            indices_t succs;
            all_successor_positions_iterative(succs, i);

            return succs;

        }

        indices_t get_predecessors(index_t j) override {

            indices_t preds;
            predecessors_init(preds, j);

            return preds;

        }

        positions_t get_range(index_t i1, index_t i2, index_t j1, index_t j2) override {

            positions_t pairs;
            range_init(pairs, i1, i2, j1, j2);

            return pairs;

        }

        bool contains_link(index_t i1, index_t i2, index_t j1, index_t j2) override {
            return link_in_range_init(i1, i2, j1, j2);
        }

        size_t count_links() override {

            size_t res = 0;
            for (size_t i = 0; i < leaves_.size(); i++) {
                res += leaves_[i];
            }

            return res;

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
            return get_range(0, n_prime_ - 1, 0, n_prime_ - 1);
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
            return link_in_range_init(i1, i2, j1, j2);
        }

        size_t count_values() override {
            return count_links();
        }


        size_t size_in_bytes() const override {
            return sizeof(BasicK2Tree<B, R, I>) // BasicK2Tree itself
                + tree_interior_.size_in_bytes() - sizeof(bit_vector_t) // tree_interior_ (do not count size of bit_vector_t object twice)
                + leaves_.size_in_bytes() - sizeof(bit_vector_t) // leaves_ (do not count size of bit_vector_t object twice)
                + ranks_.size_in_bytes() - sizeof(rank_t); // ranks_ (do not count size of rank_t object twice)
        }

        void print(bool all) override {

            std::cout << "### Parameters ###" << std::endl;
            std::cout << "h  = " << h_ << std::endl;
            std::cout << "k  = " << k_ << std::endl;
            std::cout << "n' = " << n_prime_ << std::endl;

            if (all) {

                std::cout << "### tree_interior ###" << std::endl;
                for (size_t i = 0; i < tree_interior_.size(); i++) std::cout << tree_interior_[i];
                std::cout << std::endl << std::endl;

                std::cout << "### leaves ###" << std::endl;
                for (size_t i = 0; i < leaves_.size(); i++) std::cout << leaves_[i] << " ";
                std::cout << std::endl << std::endl;

                std::cout << "### ranks ###" << std::endl;
                print_ranks(ranks_);
                std::cout << std::endl;

            }

        }

        void set_null(index_t i, index_t j, bool deep) override {

            if (deep) {
                set_null_iterative_deep(i, j);
            } else {
                set_null_iterative_shallow(i, j);
            }

        }

        void set_null_deep(index_t i, index_t j) override {
            set_null_iterative_deep(i, j);
        }

        // note: can "invalidate" the data structure (contains_link() probably won't work correctly afterwards)
        void set_null_shallow(size_t i, size_t j) override {
            set_null_iterative_shallow(i, j);
        }

        void set_null_column(index_t j, bool deep) override {

            if (deep) {
                set_null_column_iterative_deep(j);
            } else {
                set_null_column_iterative_shallow(j);
            }

        }

        void set_null_column_deep(index_t j) override {
            set_null_column_iterative_deep(j);
        }

        // note: can "invalidate" the data structure (contains_link() probably won't work correctly afterwards)
        void set_null_column_shallow(size_t j) override {
            set_null_column_iterative_shallow(j);
        }

        index_t get_first_successor(index_t i) override {
            return first_successor_position_iterative(i);
        }



    private:
        // representation of all but the last levels of the K2Tree (internal structure)
        bit_vector_t tree_interior_;

        // representation of the last level of the K2Tree (actual values of the relation)
        bit_vector_t leaves_;

        // rank data structure for navigation in tree_interior_;
        // when set_null(_column) is used in the deep version, i.e. tree_interior_ is also modified,
        // the rank structure needs to create a copy of tree_interior_ (indicated by perm_rank_ being true)
        rank_t ranks_;
        bool perm_rank_;

        size_t h_; // height of the K2Tree
        size_t k_; // arity of the K2Tree
        size_t n_prime_; // edge length of the represented relation matrix

        value_t null_ = false; // null value

        typedef RunEncodingRank<std::vector<bool>, size_t> tmp_ranks_t;


        /* helper methods for inplace construction from single list of pairs */

        size_t compute_key(const pos_t& pair, const Subproblem<size_t>& sp, size_t width) {
            return ((pair.row - sp.first_row) / width) * k_ + (pair.col - sp.first_col) / width;
        }

        void counting_sort(positions_t& pairs, std::vector<std::pair<size_t, size_t>>& intervals, const Subproblem<size_t>& sp, size_t width, size_t sup) {

            std::vector<size_t> counts(sup);

            // determine key frequencies
            for (size_t i = sp.left; i < sp.right; i++) {
                counts[compute_key(pairs[i], sp, width)]++;
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

                tmp_pairs[counts[compute_key(pairs[i], sp, width)]] = pairs[i];
                counts[compute_key(pairs[i], sp, width)]++;

            }

            for (size_t i = sp.left; i < sp.right; i++) {
                pairs[i] = tmp_pairs[i - sp.left];
            }

        }

        void build_from_pairs_inplace(positions_t& pairs) {// 3.3.5

            std::queue<Subproblem<size_t>> queue;
            std::vector<bool> tmp_interior, tmp_leaves;

            queue.emplace(0, n_prime_ - 1, 0, n_prime_ - 1, 0, pairs.size());

            while (!queue.empty()) {

                Subproblem<size_t> sp = queue.front();
                queue.pop();

                size_t sp_size = sp.last_row - sp.first_row + 1;

                if (sp_size > k_) {

                    std::vector<std::pair<size_t, size_t>> intervals(k_ * k_);
                    counting_sort(pairs, intervals, sp, sp_size / k_, k_ * k_);

                    for (size_t i = 0; i < k_ * k_; i++) {

                        if (intervals[i].first < intervals[i].second) {

                            tmp_interior.push_back(true);
                            queue.emplace(
                                    sp.first_row + (i / k_) * (sp_size / k_),
                                    sp.first_row + (i / k_ + 1) * (sp_size / k_) - 1,
                                    sp.first_col + (i % k_) * (sp_size / k_),
                                    sp.first_col + (i % k_ + 1) * (sp_size / k_) - 1,
                                    sp.left + intervals[i].first,
                                    sp.left + intervals[i].second
                            );

                        } else {
                            tmp_interior.push_back(false);
                        }

                    }

                } else {

                    std::vector<bool> leaf_values(k_ * k_, false);

                    for (size_t i = sp.left; i < sp.right; i++) {
                        leaf_values[(pairs[i].row - sp.first_row) * k_ + (pairs[i].col - sp.first_col)] = true;
                    }

                    tmp_leaves.insert(tmp_leaves.end(), leaf_values.begin(), leaf_values.end());

                }

            }

            leaves_ = bit_vector_t(tmp_leaves.size());
            leaves_.set(0, tmp_leaves);
            tmp_leaves.clear();
            tmp_leaves.shrink_to_fit();

            tree_interior_ = bit_vector_t(tmp_interior);

        }


        /* are_related() */

        bool check_link_init(size_t p, size_t q) {
            return !(leaves_.empty()) && check_link(n_prime_ / k_, p % (n_prime_ / k_), q % (n_prime_ / k_), (p / (n_prime_ / k_)) * k_ + q / (n_prime_ / k_));
        }

        bool check_link(size_t n, size_t p, size_t q, size_t z) {

            if (z >= tree_interior_.size()) {
                return leaves_[z - tree_interior_.size()];
            } else {
                return tree_interior_[z] && check_link(n / k_, p % (n / k_), q % (n / k_), ranks_.rank(z + 1) * k_ * k_ + (p / (n / k_)) * k_ + q / (n / k_));
            }

        }

        /* get_successors() */

        void all_successor_positions_iterative(indices_t& succs, size_t p) {

            if (leaves_.empty()) return;

            size_t len_interior = tree_interior_.size();

            if (len_interior == 0) {

                size_t offset = p * n_prime_;
                for (size_t i = 0; i < n_prime_; i++) {
                    if (leaves_[offset + i]) {
                        succs.push_back(i);
                    }
                }

            } else {

                std::queue<SubrowInfo<size_t>> queue, next_level_queue;

                size_t n = n_prime_/ k_;
                size_t rel_p = p;
                for (size_t j = 0, dq = 0, z = k_ * (rel_p / n); j < k_; j++, dq += n, z++) {
                    queue.emplace(dq, z);
                }


                rel_p %= n;
                n /= k_;
                for (; n > 1; rel_p %= n, n /= k_) {

                    while (!queue.empty()) {

                        auto& cur = queue.front();

                        if (tree_interior_[cur.z]) {

                            size_t y = ranks_.rank(cur.z + 1) * k_ * k_ + k_ * (rel_p / n);

                            for (size_t j = 0, new_dq = cur.dq; j < k_; j++, new_dq += n, y++) {
                                next_level_queue.emplace(new_dq, y);
                            }

                        }

                        queue.pop();

                    }

                    queue.swap(next_level_queue);

                }


                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (tree_interior_[cur.z]) {

                        size_t y = ranks_.rank(cur.z + 1) * k_ * k_ + k_ * (rel_p / n) - len_interior;

                        for (size_t j = 0, new_dq = cur.dq; j < k_; j++, new_dq += n, y++) {
                            if (leaves_[y]) {
                                succs.push_back(new_dq);
                            }
                        }

                    }

                    queue.pop();

                }

            }

        }

        /* get_first_successor() */

        size_t first_successor_position_iterative(size_t p) {

            if (leaves_.empty()) return n_prime_;

            if (tree_interior_.size() == 0) {

                size_t offset = p * n_prime_;
                for (size_t i = 0; i < n_prime_; i++) {
                    if (leaves_[offset + i]) {
                        return i;
                    }
                }

            } else {

                std::stack<ExtendedSubrowInfo<size_t>> stack;
                stack.emplace(n_prime_ / k_, n_prime_ / k_, p % (n_prime_ / k_), 0, k_ * (p / (n_prime_ / k_)), 0);

                while (!stack.empty()) {

                    auto& cur = stack.top();

                    if (cur.j == k_) {
                        stack.pop();
                    } else {

                        if (cur.z >= tree_interior_.size()) {

                            if (leaves_[cur.z - tree_interior_.size()]) {
                                return cur.dq;
                            }

                        } else {

                            if (tree_interior_[cur.z]) {
                                stack.emplace(cur.nr / k_, cur.nc / k_, cur.p % (cur.nr / k_), cur.dq, ranks_.rank(cur.z + 1) * k_ * k_ + k_ * (cur.p / (cur.nr / k_)), 0);
                            }

                        }

                        cur.dq += cur.nc;
                        cur.z++;
                        cur.j++;

                    }

                }
            }

            return n_prime_;

        }

        /* get_predecessors() */

        void predecessors_init(indices_t& preds, size_t q) {

            if (!leaves_.empty()) {

                size_t y = q / (n_prime_ / k_);

                for (size_t i = 0; i < k_; i++) {
                    predecessors(preds, n_prime_ / k_, q % (n_prime_ / k_), (n_prime_ / k_) * i, y + i * k_);
                }

            }

        }

        void predecessors(indices_t& preds, size_t n, size_t q, size_t p, size_t z) {

            if (z >= tree_interior_.size()) {

                if (leaves_[z - tree_interior_.size()]) {
                    preds.push_back(p);
                }

            } else {

                if (tree_interior_[z]) {

                    size_t y = ranks_.rank(z + 1) * k_ * k_ + q / (n / k_);

                    for (size_t i = 0; i < k_; i++) {
                        predecessors(preds, n / k_, q % (n / k_), p + (n / k_) * i, y + i * k_);
                    }

                }

            }

        }

        /* get_range() */

        void range_init(positions_t& pairs, size_t p1, size_t p2, size_t q1, size_t q2) {

            if (!leaves_.empty()) {

                for (size_t i = p1 / (n_prime_ / k_); i <= p2 / (n_prime_ / k_); i++) {

                    size_t p1_prime = (i == p1 / (n_prime_ / k_)) * (p1 % (n_prime_ / k_));
                    size_t p2_prime = (i == p2 / (n_prime_ / k_)) ? p2 % (n_prime_ / k_) : (n_prime_ / k_) - 1;

                    for (size_t j = q1 / (n_prime_ / k_); j <= q2 / (n_prime_ / k_); j++) {
                        range(
                                pairs,
                                n_prime_ / k_,
                                p1_prime,
                                p2_prime,
                                (j == q1 / (n_prime_ / k_)) * (q1 % (n_prime_ / k_)),
                                (j == q2 / (n_prime_ / k_)) ? q2 % (n_prime_ / k_) : (n_prime_ / k_) - 1,
                                (n_prime_ / k_) * i,
                                (n_prime_ / k_) * j,
                                k_ * i + j
                        );
                    }

                }

            }

        }

        void range(positions_t& pairs, size_t n, size_t p1, size_t p2, size_t q1, size_t q2, size_t dp, size_t dq, size_t z) {

            if (z >= tree_interior_.size()) {

                if (leaves_[z - tree_interior_.size()]) {
                    pairs.emplace_back(dp, dq);
                }

            } else {

                if (tree_interior_[z]) {

                    size_t y = ranks_.rank(z + 1) * k_ * k_;

                    for (size_t i = p1 / (n / k_); i <= p2 / (n / k_); i++) {

                        size_t p1_prime = (i == p1 / (n / k_)) * (p1 % (n / k_));
                        size_t p2_prime = (i == p2 / (n / k_)) ? p2 % (n / k_) : n / k_ - 1;

                        for (size_t j = q1 / (n / k_); j <= q2 / (n / k_); j++) {
                            range(
                                    pairs,
                                    n / k_,
                                    p1_prime,
                                    p2_prime,
                                    (j == q1 / (n / k_)) * (q1 % (n / k_)),
                                    (j == q2 / (n / k_)) ? q2 % (n / k_) : n / k_ - 1,
                                    dp + (n / k_) * i,
                                    dq + (n / k_) * j,
                                    y + k_ * i + j
                            );
                        }

                    }

                }

            }

        }

        /* link_in_range() */

        bool link_in_range_init(size_t p1, size_t p2, size_t q1, size_t q2) {

            if (!leaves_.empty()) {

                // dividing by k_ (as stated in the paper) in not correct,
                // because it does not use the size of the currently considered submatrix but of its submatrices
                if ((p1 == 0) && (q1 == 0) && (p2 == (n_prime_ /*/ k_*/ - 1)) && (q2 == (n_prime_ /*/ k_*/ - 1))) {
                    return true;
                }

                for (size_t i = p1 / (n_prime_ / k_); i <= p2 / (n_prime_ / k_); i++) {

                    size_t p1_prime = (i == p1 / (n_prime_ / k_)) * (p1 % (n_prime_ / k_));
                    size_t p2_prime = (i == p2 / (n_prime_ / k_)) ? (p2 % (n_prime_ / k_)) : n_prime_ / k_ - 1;

                    for (size_t j = q1 / (n_prime_ / k_); j <= q2 / (n_prime_ / k_); j++) {

                        if (link_in_range(n_prime_ / k_, p1_prime, p2_prime, (j == q1 / (n_prime_ / k_)) * (q1 % (n_prime_ / k_)), (j == q2 / (n_prime_ / k_)) ? q2 % (n_prime_ / k_) : n_prime_ / k_ - 1, k_ * i + j)) {
                            return true;
                        }

                    }

                }

            }

            return false;

        }

        bool link_in_range(size_t n, size_t p1, size_t p2, size_t q1, size_t q2, size_t z) {

            if (z >= tree_interior_.size()) {

                return leaves_[z - tree_interior_.size()];

            } else {

                if (tree_interior_[z]) {

                    // dividing by k_ (as stated in the paper) in not correct,
                    // because it does not use the size of the currently considered submatrix but of its submatrices
                    if ((p1 == 0) && (q1 == 0) && (p2 == (n /*/ k_*/ - 1)) && (q2 == (n /*/ k_*/ - 1))) {
                        return true;
                    }

                    size_t y = ranks_.rank(z + 1) * k_ * k_;

                    for (size_t i = p1 / (n / k_); i <= p2 / (n / k_); i++) {

                        size_t p1_prime = (i == p1 / (n / k_)) * (p1 % (n / k_));
                        size_t p2_prime = (i == p2 / (n / k_)) ? (p2 % (n / k_)) : n / k_ - 1;

                        for (size_t j = q1 / (n / k_); j <= q2 / (n / k_); j++) {

                            if (link_in_range(n / k_, p1_prime, p2_prime, (j == q1 / (n / k_)) * (q1 % (n / k_)), (j == q2 / (n / k_)) ? q2 % (n / k_) : n / k_ - 1, y + k_ * i + j)) {
                                return true;
                            }

                        }

                    }

                }

                return false;

            }

        }


        /* set_null() */

        void set_null_iterative_shallow(size_t i, size_t j) {

            if (leaves_.empty()) return;

            size_t z = (i / (n_prime_ / k_)) * k_ + j / (n_prime_ / k_);
            size_t p = i % (n_prime_ / k_);
            size_t q = j % (n_prime_ / k_);
            size_t n = n_prime_ / k_;

            for (; z < tree_interior_.size() && tree_interior_[z]; ) {

                z = ranks_.rank(z + 1) * k_ * k_ + (p / (n / k_)) * k_ + q / (n / k_);
                p = p % (n / k_);
                q = q % (n / k_);
                n = n / k_;

            }

            if (z >= tree_interior_.size()) {
                leaves_.set(z - tree_interior_.size(), false);
            }

        }

        void set_null_iterative_deep(size_t i, size_t j) {

            if (leaves_.empty()) return;

            std::stack<std::pair<size_t, size_t>> pos_stack;

            size_t z = (i / (n_prime_ / k_)) * k_ + j / (n_prime_ / k_);
            size_t p = i % (n_prime_ / k_);
            size_t q = j % (n_prime_ / k_);
            size_t n = n_prime_ / k_;

            pos_stack.emplace(0, z);

            for (; z < tree_interior_.size() && tree_interior_[z]; ) {

                pos_stack.emplace(ranks_.rank(z + 1) * k_ * k_, (p / (n / k_)) * k_ + q / (n / k_));

                z = pos_stack.top().first + pos_stack.top().second;
                p = p % (n / k_);
                q = q % (n / k_);
                n = n / k_;

            }

            if (z >= tree_interior_.size() && leaves_[z - tree_interior_.size()]) {

                size_t bs = pos_stack.top().first;
                size_t os = pos_stack.top().second;
                pos_stack.pop();
                leaves_.set((bs + os) - tree_interior_.size(), false);

                for (bool cont = leaves_.is_zero(bs - tree_interior_.size(), k_ * k_);
                     !pos_stack.empty() && cont;
                     cont = tree_interior_.is_zero(bs, k_ * k_)) {

                    bs = pos_stack.top().first;
                    os = pos_stack.top().second;
                    pos_stack.pop();
                    tree_interior_.set(bs + os, false);

                }

            }

        }

        void set_null_column_iterative_shallow(size_t q) {

            if (leaves_.empty()) return;

            size_t len_interior = tree_interior_.size();

            if (len_interior == 0) {

                for (size_t i = 0; i < n_prime_; i++) {
                    if (leaves_[i * n_prime_ + q] != null_) {
                        leaves_.set(i * n_prime_ + q, false);
                    }
                }

            } else {

                std::queue<SubcolInfo<size_t>> queue, next_level_queue;

                size_t n = n_prime_/ k_;
                size_t rel_q = q;
                for (size_t i = 0, dp = 0, z = q / (n_prime_ / k_); i < k_; i++, dp += n, z += k_) {
                    queue.emplace(dp, z);
                }

                rel_q %= n;
                n /= k_;
                for (; n > 1; rel_q %= n, n /= k_) {

                    while (!queue.empty()) {

                        auto& cur = queue.front();

                        if (tree_interior_[cur.z]) {

                            size_t y = ranks_.rank(cur.z + 1) * k_ * k_ + rel_q / n;

                            for (size_t j = 0, new_dp = cur.dp; j < k_; j++, new_dp += n, y += k_) {
                                next_level_queue.emplace(new_dp, y);
                            }

                        }

                        queue.pop();

                    }

                    queue.swap(next_level_queue);

                }

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (tree_interior_[cur.z]) {

                        size_t y = ranks_.rank(cur.z + 1) * k_ * k_ + rel_q / n - len_interior;

                        for (size_t j = 0, new_dp = cur.dp; j < k_; j++, new_dp += n, y += k_) {
                            if (leaves_[y] != null_) {
                                leaves_.set(y, null_);
                            }
                        }

                    }

                    queue.pop();

                }

            }

        }

        void set_null_column_iterative_deep(size_t q) {

            if (leaves_.empty()) return;

            size_t len_interior = tree_interior_.size();

            if (len_interior == 0) {

                for (size_t i = 0; i < n_prime_; i++) {
                    if (leaves_[i * n_prime_ + q] != null_) {
                        leaves_.set(i * n_prime_ + q, false);
                    }
                }

            } else {

                std::stack<SubcolInfoDeep<size_t>> stack;
                bool sideways = false;
                bool all_zero = true;

                stack.emplace(0, q, 0, q / (n_prime_ / k_), n_prime_ / k_, n_prime_ / k_, 0, 0, true);

                while (!stack.empty()) {

                    auto& cur = stack.top();

                    if (cur.nr == 1) { // leaves

                        size_t y = cur.z - len_interior;
                        all_zero = false;

                        for (size_t j = 0; j < k_; j++, y += k_) {

                            all_zero |= leaves_[y] != null_;
                            leaves_.set(y, null_);

                        }

                        // initially at least one non-null value in each leaf, call is_all() only if at least one non-null value was eliminated
                        all_zero = all_zero && leaves_.is_zero(cur.node - len_interior, k_ * k_);

                        stack.pop();
                        sideways = true;

                    } else { // interior

                        if (sideways) {

                            tree_interior_.set(cur.z, !all_zero);
                            cur.zero &= all_zero;
                            cur.dp += cur.nr;
                            cur.z += k_;
                            cur.i += 1;

                        }

                        while (cur.i < k_ && !tree_interior_[cur.z]) {

                            cur.dp += cur.nr;
                            cur.z += k_;
                            cur.i += 1;

                        }

                        if (cur.i < k_) {

                            stack.emplace(cur.dp, cur.dq % cur.nc, ranks_.rank(cur.z + 1) * k_ * k_,
                                          (cur.dq % cur.nc) / (cur.nc / k_), cur.nr / k_, cur.nc / k_, 0, 0, true);
                            sideways = false;

                        } else {

                            all_zero = cur.zero && tree_interior_.is_zero(cur.node, k_ * k_);
                            stack.pop();
                            sideways = true;

                        }

                    }

                }

            }

        }

    };

}

#endif //K2TREES_BASICK2TREE_HPP
