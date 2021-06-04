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

#ifndef GEFAST_SPACEPARTITIONING_HPP
#define GEFAST_SPACEPARTITIONING_HPP

#include <fstream>
#include <numeric>
#include <sstream>

#include "Base.hpp"

namespace GeFaST {

    // available norm tags used to adjust the behaviour of k-d trees
    struct manhattan_norm_tag { }; // Manhattan norm / distance
    struct euclidean_norm_tag { }; // Euclidean norm / distance

    /*
     * Data structure for dividing a(n Euclidean) space into non-overlapping region.
     * Allows efficient geometric queries (e.g. finding a nearest object or all within a radius) on stored objects.
     */
    struct SpacePartitioningStructure {

        virtual ~SpacePartitioningStructure() = default;

        /* === nearest-neighbour methods */

        /*
         * Find the nearest neighbour to the input element with index 'query'.
         *
         * Prevents a self match.
         */
        virtual numSeqs_t nearest(numSeqs_t query) const = 0;

        /*
         * Find the nearest neighbour to the given coordinates.
         */
        virtual numSeqs_t nearest_value(const feat_t* query) const = 0;

        /*
         * Find the nearest neighbour to the input element with index 'query'.
         * Return the index of that neighbour and its distance.
         *
         * Prevents a self match.
         */
        virtual AuxiliaryData::Partner nearest_partner(numSeqs_t query) const = 0;

        /*
         * Find the nearest neighbour to the given coordinates.
         * Return the index of that neighbour and its distance.
         */
        virtual AuxiliaryData::Partner nearest_partner_value(const feat_t* query) const = 0;


        /* === range-search methods === */

        /*
         * Find all elements that have a distance at most 'rad'
         * from the input element with index 'query'.
         *
         * Prevents a self match.
         */
        virtual std::vector<numSeqs_t> range(numSeqs_t query, dist_t rad) const = 0;

        /*
         * Find all elements that have a distance at most 'rad'
         * to the given coordinates.
         */
        virtual std::vector<numSeqs_t> range_value(const feat_t* query, dist_t rad) const = 0;

        /*
         * Find all elements that have a distance at most 'rad' from the input element with index 'query'.
         * Return the indices of these neighbour and their distances.
         *
         * Prevents a self match.
         */
        virtual std::vector<AuxiliaryData::Partner> range_partner(numSeqs_t query, dist_t rad) const = 0;

        /*
         * Find all elements that have a distance at most 'rad'
         * to the given coordinates.
         */
        virtual std::vector<AuxiliaryData::Partner> range_partner_value(const feat_t* query, dist_t rad) const = 0;

    };




    /*
     * Implementation of k-d trees allowing for nearest neighbour and range search.
     * The k-d tree is built directly from the given coordinates (features).
     *
     * Currently, a balanced k-d tree is constructed by determining the actual medians during tree construction.
     *
     * The search methods return integer indices as follows:
     * If an integer i is returned, the element whose coordinates were given to the constructor at index i
     * is the nearest neighbour resp. in the range.
     *
     * Template parameters:
     *  - N = tag of used norm (see above)
     *
     *  TODO? approximate median-selection method (faster construction vs potentially slower queries (less balanced))
     *  TODO? lazy delete (to mark swarmed amplicons, avoids removing them afterwards in AuxiliaryData)
     *  TODO? versions of query methods with additional parameters to check abundance (if desired) and swarmed status
     *     (avoids removing the corresponding amplicons afterwards in AuxiliaryData)
     *
     * Based on:
     * "Multidimensional Binary Search Trees Used for Associative Searching"
     * (Bentley, 1975, https://doi.org/10.1145/361002.361007)
     * with code adapted and extended from https://rosettacode.org/wiki/K-d_tree.
     */
    template<typename N>
    class KdTree : public SpacePartitioningStructure {

    public:
        KdTree(size_t dim, const feat_t* coordinates, numSeqs_t num, bool own = false) {

            dim_ = dim;
            owns_coordinates_ = own;

            num_nodes_ = num;
            num_elements_ = num;

            nodes_ = new Node[num_nodes_ + 1]; // + 1 to include "error position" and avoid extra cases
            coordinates_ = coordinates;

            // create nodes of k-d tree (not yet connected)
            for (numSeqs_t i = 0; i < num_nodes_; i++) {
                nodes_[i].id = i;
            }
            nodes_[num_nodes_].id = nodes_[num_nodes_].left = nodes_[num_nodes_].right = num_nodes_;

            // build actual k-d tree from features
            pos_root_ = make_tree(0, num_nodes_, 0);

        }

        KdTree(size_t dim, const feat_t* coordinates, numSeqs_t num, const std::vector<numSeqs_t>& selection, bool own = false) {

            dim_ = dim;
            owns_coordinates_ = own;

            num_nodes_ = selection.size();
            num_elements_ = num;

            nodes_ = new Node[num_nodes_ + 1]; // + 1 to include "error position" and avoid extra cases
            coordinates_ = coordinates;

            // create nodes of k-d tree (not yet connected)
            for (numSeqs_t i = 0; i < num_nodes_; i++) {
                nodes_[i].id = selection[i];
            }
            nodes_[num_nodes_].id = nodes_[num_nodes_].left = nodes_[num_nodes_].right = num_nodes_;

            // build actual k-d tree from features
            pos_root_ = make_tree(0, num_nodes_, 0);

        }

        explicit KdTree(const FeatureAmpliconCollection& ac)
            : KdTree(ac.num_features(), ac.all_features(), ac.size(), false) {
            // nothing else to do
        }

        KdTree(const FeatureAmpliconCollection& ac, const std::vector<numSeqs_t>& selection)
            : KdTree(ac.num_features(), ac.all_features(), ac.size(), selection, false) {
            // nothing else to do
        }

        ~KdTree() {

            if (owns_coordinates_) delete[] coordinates_;
            delete[] nodes_;

        }

        KdTree* clone() const { // deep-copy clone method
            return new KdTree(*this);
        }

        KdTree(const KdTree& other) { // copy constructor

            dim_ = other.dim_;

            num_nodes_ = other.num_nodes_;
            num_elements_ = other.num_elements_;

            nodes_ = new Node[num_nodes_ + 1]; // + 1 to include "error position" and avoid extra cases
            std::copy(other.nodes_, other.nodes_ + num_nodes_ + 1, nodes_);

            if (other.owns_coordinates_) {

                auto tmp = new feat_t[num_elements_ * dim_];
                std::copy(other.coordinates_, other.coordinates_ + num_elements_ * dim_, tmp);
                coordinates_ = tmp;

            } else {
                coordinates_ = other.coordinates_;
            }

            owns_coordinates_ = other.owns_coordinates_;
            pos_root_ = other.pos_root_;

        }

        KdTree(KdTree&& other) noexcept { // move constructor

            dim_ = other.dim_;

            num_nodes_ = other.num_nodes_;
            num_elements_ = other.num_elements_;

            nodes_ = other.nodes_; other.nodes_ = nullptr;
            coordinates_ = other.coordinates_; other.coordinates_ = nullptr;

            owns_coordinates_ = other.owns_coordinates_;
            pos_root_ = other.pos_root_;

        }

        KdTree& operator=(const KdTree& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            if (owns_coordinates_) delete[] coordinates_;
            delete[] nodes_;

            // copy new resources
            dim_ = other.dim_;

            num_nodes_ = other.num_nodes_;
            num_elements_ = other.num_elements_;

            nodes_ = new Node[num_nodes_ + 1]; // + 1 to include "error position" and avoid extra cases
            std::copy(other.nodes_, other.nodes_ + num_nodes_ + 1, nodes_);

            if (other.owns_coordinates_) {

                auto tmp = new feat_t[num_elements_ * dim_];
                std::copy(other.coordinates_, other.coordinates_ + num_elements_ * dim_, tmp);
                coordinates_ = tmp;

            } else {
                coordinates_ = other.coordinates_;
            }

            owns_coordinates_ = other.owns_coordinates_;
            pos_root_ = other.pos_root_;

            return *this;

        }

        KdTree& operator=(KdTree&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            if (owns_coordinates_) delete[] coordinates_;
            delete[] nodes_;

            // move new resources
            dim_ = other.dim_;

            num_nodes_ = other.num_nodes_;
            num_elements_ = other.num_elements_;

            nodes_ = other.nodes_; other.nodes_ = nullptr;
            coordinates_ = other.coordinates_; other.coordinates_ = nullptr;

            owns_coordinates_ = other.owns_coordinates_;
            pos_root_ = other.pos_root_;

            return *this;

        }


        /* === nearest-neighbour methods */

        /*
         * Find the nearest neighbour to the input element with index 'query'.
         *
         * Prevents a self match.
         */
        numSeqs_t nearest(numSeqs_t query) const override {

            numSeqs_t best = num_nodes_;
            dist_t best_dist = 0;
            auto coords = coordinates_ + query * dim_;

            nearest(pos_root_, coords, 0, &best, &best_dist, query);

            return nodes_[best].id;

        }

        /*
         * Find the nearest neighbour to the given coordinates.
         */
        numSeqs_t nearest_value(const feat_t* query) const override {

            numSeqs_t best = num_nodes_;
            dist_t best_dist = 0;

            nearest(pos_root_, query, 0, &best, &best_dist);

            return nodes_[best].id;

        }

        /*
         * Find the nearest neighbour to the input element with index 'query'.
         * Return the index of that neighbour and its distance.
         *
         * Prevents a self match.
         */
        AuxiliaryData::Partner nearest_partner(numSeqs_t query) const override {

            numSeqs_t best = num_nodes_;
            dist_t best_dist = 0;
            auto coords = coordinates_ + query * dim_;

            nearest(pos_root_, coords, 0, &best, &best_dist, query);

            return AuxiliaryData::Partner(nodes_[best].id, get_root(best_dist, norm_tag_));

        }

        /*
         * Find the nearest neighbour to the given coordinates.
         * Return the index of that neighbour and its distance.
         */
        AuxiliaryData::Partner nearest_partner_value(const feat_t* query) const override {

            numSeqs_t best = num_nodes_;
            dist_t best_dist = 0;

            nearest(pos_root_, query, 0, &best, &best_dist);

            return AuxiliaryData::Partner(nodes_[best].id, get_root(best_dist, norm_tag_));

        }


        /* === range-search methods === */

        /*
         * Find all elements that have a distance at most 'rad'
         * from the input element with index 'query'.
         *
         * Prevents a self match.
         */
        std::vector<numSeqs_t> range(numSeqs_t query, dist_t rad) const override {

            auto coords = coordinates_ + query * dim_;
            std::vector<numSeqs_t> rng;

            range(pos_root_, coords, rad, get_term(rad, norm_tag_), 0, rng, query);

            return rng;

        }

        /*
         * Find all elements that have a distance at most 'rad'
         * to the given coordinates.
         */
        std::vector<numSeqs_t> range_value(const feat_t* query, dist_t rad) const override {

            std::vector<numSeqs_t> rng;

            range(pos_root_, query, rad, get_term(rad, norm_tag_), 0, rng);

            return rng;

        }

        /*
         * Find all elements that have a distance at most 'rad' from the input element with index 'query'.
         * Return the indices of these neighbour and their distances.
         *
         * Prevents a self match.
         */
        std::vector<AuxiliaryData::Partner> range_partner(numSeqs_t query, dist_t rad) const override {

            auto coords = coordinates_ + query * dim_;
            std::vector<AuxiliaryData::Partner> rng;

            range(pos_root_, coords, rad, get_term(rad, norm_tag_), 0, rng, query);

            return rng;

        }

        /*
         * Find all elements that have a distance at most 'rad'
         * to the given coordinates.
         */
        std::vector<AuxiliaryData::Partner> range_partner_value(const feat_t* query, dist_t rad) const override {

            std::vector<AuxiliaryData::Partner> rng;

            range(pos_root_, query, rad, get_term(rad, norm_tag_), 0, rng);

            return rng;

        }

    protected:
        /*
         * Default constructor -- not to be used alone, only for derived classes.
         */
        KdTree() {

            dim_ = 0;

            num_nodes_ = 0;
            num_elements_ = 0;

            nodes_ = nullptr;
            coordinates_ = nullptr;
            owns_coordinates_ = false;

            pos_root_ = 0;

        }

        /*
         * Return coordinate of d-th dimension of element stored at index a in nodes_.
         */
        inline feat_t get_coord(numSeqs_t a, size_t d) const {
            return coordinates_[nodes_[a].id * dim_ + d];
        }

        /*
         * Swap the two nodes at the given indices in nodes_.
         */
        void swap(numSeqs_t a, numSeqs_t b) {

            std::swap(nodes_[a].id, nodes_[b].id);
            std::swap(nodes_[a].left, nodes_[b].left);
            std::swap(nodes_[a].right, nodes_[b].right);

        }

        /*
         * Compute the distance (without taking the root) between the elements
         * represented by the nodes at the two given indices in nodes_.
         */
        dist_t distance(numSeqs_t a, numSeqs_t b) {

            dist_t sum = 0;
            for (auto d = 0; d < dim_; d++) {

                auto tmp = get_coord(a, d) - get_coord(b, d);
                sum += get_term(tmp, norm_tag_);

            }

            return sum;

        }

        /*
         * Compute the distance (without taking the root) between the node
         * at the given index in nodes_ and another position specified by its coordinates.
         */
        dist_t distance(numSeqs_t a, const feat_t* coords) const {

            dist_t sum = 0;
            for (auto d = 0; d < dim_; d++) {

                auto tmp = get_coord(a, d) - coords[d];
                sum += get_term(tmp, norm_tag_);

            }

            return sum;

        }

        /*
         * Auxiliary method of 'make_tree' for finding the position of the exact median element in the range [begin, end)
         * and partially sorting the range.
         *
         * Iterative quickselect (less before median, greater-equal from median onwards).
         */
        numSeqs_t find_median(numSeqs_t begin, numSeqs_t end, size_t d) {

            numSeqs_t mid = begin + (end - begin) / 2;

            while (true) {

                if (end == begin + 1) return begin;

                feat_t pivot = get_coord(mid, d);
                numSeqs_t store = begin;

                swap(mid, end - 1);
                for (numSeqs_t p = begin; p < end; p++) {

                    if (get_coord(p, d) < pivot) {

                        if (p != store) swap(p, store);
                        store++;

                    }

                }
                swap(store, end - 1);

                if (mid == store) return mid;

                if (mid < store) {
                    end = store;
                } else {
                    begin = store + 1;
                }

            }

        }

        /*
         * Recursively construct the k-d tree from the nodes using their coordinates.
         */
        numSeqs_t make_tree(numSeqs_t begin, numSeqs_t len, size_t d) {

            numSeqs_t n;

            if (len == 0) return num_nodes_;

            if ((n = find_median(begin, begin + len, d)) < num_nodes_) {

                d = (d + 1) % dim_;
                nodes_[n].left  = make_tree(begin, n - begin, d);
                nodes_[n].right = make_tree(n + 1, begin + len - (n + 1), d);

            }

            return n;

        }

        /*
         * Auxiliary method providing the key functionality to find the nearest neighbour of a specified position.
         *
         * Parameters:
         *  - root: root node of the current subtree to search for the nearest neighbour
         *  - query: coordinates of the position relative to which the nearest neighbour is searched
         *  - d: splitting dimension of the current node ('root')
         *  - best: index (in 'nodes_') of the currently nearest neighbour
         *  - best_dist: distance of current 'best' to 'query'
         */
        void nearest(numSeqs_t root, const feat_t* query, size_t d, numSeqs_t* best, dist_t* best_dist) const {

            if (root == num_nodes_) return; // root == num_nodes_ indicates that current subtree is empty

            dist_t dist = distance(root, query);
            dist_t dist_dim = get_coord(root, d) - query[d];

            // update when nearer neighbour found (or with the first neighbour encountered (-> second comparison))
            if ((dist < *best_dist) || (*best == num_nodes_)) {

                *best_dist = dist;
                *best = root;

            }

            // directly return with current best when it is equal (or, w.r.t. precision, practically equal)
            if (!*best_dist) return;

            // looking only at the current dimension,
            // decide whether the query is smaller than the current root (-> side of the query is the left subtree)
            // or at least as large (-> side of the query is the right subtree)
            // [fits partitioning in make_tree / find_median]
            numSeqs_t query_side, other_side;
            if (dist_dim > 0) {

                query_side = nodes_[root].left;
                other_side = nodes_[root].right;

            } else {

                query_side = nodes_[root].right;
                other_side = nodes_[root].left;

            }

            d = (d + 1) % dim_;

            nearest(query_side, query, d, best, best_dist);
            // there can be nearer neighbour on the other side of the current hyperplane
            // only if distance implied by the current distance does not already meet or exceed the currently best distance
            if (get_term(dist_dim, norm_tag_) < *best_dist) {
                nearest(other_side, query, d, best, best_dist);
            }

        }

        /*
         * Auxiliary method providing the key functionality to find the nearest neighbour of a specified position
         * and avoiding self matches by not allowing the identifier of the found node to equal 'query_id'.
         *
         * Other parameters as above.
         */
        void nearest(numSeqs_t root, const feat_t* query, size_t d, numSeqs_t* best, dist_t* best_dist, numSeqs_t query_id) const {

            if (root == num_nodes_) return; // root == num_nodes_ indicates that current subtree is empty

            dist_t dist = distance(root, query);
            dist_t dist_dim = get_coord(root, d) - query[d];

            // update when nearer neighbour found (or with the first neighbour encountered (-> second comparison))
            if (((dist < *best_dist) || (*best == num_nodes_)) && (nodes_[root].id != query_id)) {

                *best_dist = dist;
                *best = root;

            }

            // directly return with current best when it is equal (or, w.r.t. precision, practically equal)
            if (!*best_dist) return;

            // looking only at the current dimension,
            // decide whether the query is smaller than the current root (-> side of the query is the left subtree)
            // or at least as large (-> side of the query is the right subtree)
            // [fits partitioning in make_tree / find_median]
            numSeqs_t query_side, other_side;
            if (dist_dim > 0) {

                query_side = nodes_[root].left;
                other_side = nodes_[root].right;

            } else {

                query_side = nodes_[root].right;
                other_side = nodes_[root].left;

            }

            d = (d + 1) % dim_;

            nearest(query_side, query, d, best, best_dist, query_id);
            // there can be nearer neighbour on the other side of the current hyperplane
            // only if distance implied by the current distance does not already meet or exceed the currently best distance
            if (get_term(dist_dim, norm_tag_) < *best_dist) {
                nearest(other_side, query, d, best, best_dist, query_id);
            }

        }

        /*
         * Auxiliary method providing the key functionality to find all elements with a distance at most 'rad'
         * to the given coordinates.
         *
         * Parameters:
         *  - root: root node of the current subtree to search for the nearest neighbour
         *  - query: coordinates of the position relative to which the nearest neighbour is searched
         *  - rad: search radius of the range
         *  - hd_rad: 'rad' to the power p where p is the parameter of the used norm (e.g. 2 for Euclidean norm)
         *  - d: splitting dimension of the current node ('root')
         *  - rng: collection of already found elements in the search range
         *
         * 'hd_rad' is provided to avoid computing the root in the distance method
         * or the power of the radius in each recursive call.
         *
         * Used by range_value(const feat_t*, dist_t).
         */
        void range(numSeqs_t root, const feat_t* query, dist_t rad, dist_t hd_rad, size_t d, std::vector<numSeqs_t>& rng) const {

            if (root != num_nodes_) { // root == num_nodes_ indicates that current subtree is empty

                dist_t dist = distance(root, query);
                dist_t dist_dim = get_coord(root, d) - query[d];

                if (dist <= hd_rad) { // element of current tree node falls within range
                    rng.push_back(nodes_[root].id);
                }

                // looking only at the current dimension,
                // decide whether the query is smaller than the current root (-> side of the query is the left subtree)
                // or at least as large (-> side of the query is the right subtree)
                // [fits partitioning in make_tree / find_median]
                numSeqs_t query_side, other_side;
                if (dist_dim > 0) {

                    query_side = nodes_[root].left;
                    other_side = nodes_[root].right;

                } else {

                    query_side = nodes_[root].right;
                    other_side = nodes_[root].left;

                }

                d = (d + 1) % dim_;
                range(query_side, query, rad, hd_rad, d, rng);
                // there can be elements on the other side of the current hyperplane (but still within the searched range)
                // only if the hypersphere around the query intersects with or touches the hyperplane
                if (std::abs(dist_dim) <= rad) {
                    range(other_side, query, rad, hd_rad, d, rng);
                }

            }

        }

        /*
         * Auxiliary method providing the key functionality to find all elements with a distance at most 'rad' to the given
         * coordinates and avoiding self matches by not allowing the identifiers of the found node to equal 'query_id'.
         *
         * Other parameters as above.
         *
         * Used by range(numSeqs_t, dist_t).
         */
        void range(numSeqs_t root, const feat_t* query, dist_t rad, dist_t hd_rad, size_t d, std::vector<numSeqs_t>& rng, numSeqs_t query_id) const {

            if (root != num_nodes_) { // root == num_nodes_ indicates that current subtree is empty

                dist_t dist = distance(root, query);
                dist_t dist_dim = get_coord(root, d) - query[d];

                if ((dist <= hd_rad) && (nodes_[root].id != query_id)) { // element of current tree node falls within range
                    rng.push_back(nodes_[root].id);
                }

                // looking only at the current dimension,
                // decide whether the query is smaller than the current root (-> side of the query is the left subtree)
                // or at least as large (-> side of the query is the right subtree)
                // [fits partitioning in make_tree / find_median]
                numSeqs_t query_side, other_side;
                if (dist_dim > 0) {

                    query_side = nodes_[root].left;
                    other_side = nodes_[root].right;

                } else {

                    query_side = nodes_[root].right;
                    other_side = nodes_[root].left;

                }

                d = (d + 1) % dim_;
                range(query_side, query, rad, hd_rad, d, rng, query_id);
                // there can be elements on the other side of the current hyperplane (but still within the searched range)
                // only if the hypersphere around the query intersects with or touches the hyperplane
                if (std::abs(dist_dim) <= rad) {
                    range(other_side, query, rad, hd_rad, d, rng, query_id);
                }

            }

        }

        /*
         * Auxiliary method providing the key functionality to find all elements with a distance at most 'rad'
         * to the given coordinates. Return the indices of these neighbour and their distances.
         *
         * Other parameters as above.
         *
         * Used by range_partner_value(const feat_t*, dist_t).
         */
        void range(numSeqs_t root, const feat_t* query, dist_t rad, dist_t hd_rad, size_t d, std::vector<AuxiliaryData::Partner>& rng) const {

            if (root != num_nodes_) { // root == num_nodes_ indicates that current subtree is empty

                dist_t dist = distance(root, query);
                dist_t dist_dim = get_coord(root, d) - query[d];

                if (dist <= hd_rad) { // element of current tree node falls within range
                    rng.emplace_back(nodes_[root].id, get_root(dist, norm_tag_));
                }

                // looking only at the current dimension,
                // decide whether the query is smaller than the current root (-> side of the query is the left subtree)
                // or at least as large (-> side of the query is the right subtree)
                // [fits partitioning in make_tree / find_median]
                numSeqs_t query_side, other_side;
                if (dist_dim > 0) {

                    query_side = nodes_[root].left;
                    other_side = nodes_[root].right;

                } else {

                    query_side = nodes_[root].right;
                    other_side = nodes_[root].left;

                }

                d = (d + 1) % dim_;
                range(query_side, query, rad, hd_rad, d, rng);
                // there can be elements on the other side of the current hyperplane (but still within the searched range)
                // only if the hypersphere around the query intersects with or touches the hyperplane
                if (std::abs(dist_dim) <= rad) {
                    range(other_side, query, rad, hd_rad, d, rng);
                }

            }

        }

        /*
         * Auxiliary method providing the key functionality to find all elements with a distance at most 'rad' to the given
         * coordinates and avoiding self matches by not allowing the identifiers of the found node to equal 'query_id'.
         * Return the indices of these neighbour and their distances.
         *
         * Other parameters as above.
         *
         * Used by range_partner(numSeqs_t, dist_t).
         */
        void range(numSeqs_t root, const feat_t* query, dist_t rad, dist_t hd_rad, size_t d, std::vector<AuxiliaryData::Partner>& rng, numSeqs_t query_id) const {

            if (root != num_nodes_) { // root == num_nodes_ indicates that current subtree is empty

                dist_t dist = distance(root, query);
                dist_t dist_dim = get_coord(root, d) - query[d];

                if (dist <= hd_rad && nodes_[root].id != query_id) { // element of current tree node falls within range
                    rng.emplace_back(nodes_[root].id, get_root(dist, norm_tag_));
                }

                // looking only at the current dimension,
                // decide whether the query is smaller than the current root (-> side of the query is the left subtree)
                // or at least as large (-> side of the query is the right subtree)
                // [fits partitioning in make_tree / find_median]
                numSeqs_t query_side, other_side;
                if (dist_dim > 0) {

                    query_side = nodes_[root].left;
                    other_side = nodes_[root].right;

                } else {

                    query_side = nodes_[root].right;
                    other_side = nodes_[root].left;

                }

                d = (d + 1) % dim_;
                range(query_side, query, rad, hd_rad, d, rng, query_id);
                // there can be elements on the other side of the current hyperplane (but still within the searched range)
                // only if the hypersphere around the query intersects with or touches the hyperplane
                if (std::abs(dist_dim) <= rad) {
                    range(other_side, query, rad, hd_rad, d, rng, query_id);
                }

            }

        }


        /*
         * Auxiliary class for representation of a node in a k-d tree.
         * Member 'id' refers to the input index of the coordinates represented by this node.
         * Members 'left' and 'right' point to the root of the left resp. right subtree
         * by describing their index in the 'nodes_' array.
         */
        struct Node {

            numSeqs_t id = 0;
            numSeqs_t left = 0;
            numSeqs_t right = 0;

        };

        Node* nodes_; // representation of k-d tree as list of "pointered" nodes
        numSeqs_t num_nodes_; // number of nodes in k-d tree
        numSeqs_t pos_root_; // index of tree root in nodes_

        const feat_t* coordinates_; // coordinates (features) of the elements represented in the k-d tree
        numSeqs_t num_elements_; // number of "feature vectors" in coordinates_
        bool owns_coordinates_; // flag indicating whether the coordinates_ array is owned by the k-d tree
        size_t dim_; // number of dimensions (features) per element


        N norm_tag_; // tag of used norm (influences how the dimension-wise terms and final distance are computed)

        inline dist_t get_term(dist_t val, manhattan_norm_tag t) const {
            return std::abs(val);
        }
        inline dist_t get_root(dist_t val, manhattan_norm_tag t) const {
            return val;
        }

        inline dist_t get_term(dist_t val, euclidean_norm_tag t) const {
            return val * val;
        }
        inline dist_t get_root(dist_t val, euclidean_norm_tag t) const {
            return std::sqrt(val);
        }

    };

}

#endif //GEFAST_SPACEPARTITIONING_HPP
