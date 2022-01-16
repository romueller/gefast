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

#ifndef GEFAST_FLEXIBLERELATIONS_HPP
#define GEFAST_FLEXIBLERELATIONS_HPP

#include <algorithm>

#include "../AuxiliaryData.hpp"
#include "../Base.hpp"
#include "../Utility.hpp"
#include "../modes/SpaceLevenshteinMode.hpp"
#include "k2trees/K2Tree.hpp"
#include "Basics.hpp"

namespace GeFaST {

    /*
     * Helper data structure mapping an arbitrary (non-consecutive) ascending sequence
     * of n unique (non-negative) integers onto [0:n-1].
     */
    class RankedAscendingLabels {

    public:
        RankedAscendingLabels();

        /*
         * Prepare the data structure for storing the provided number of integers.
         */
        RankedAscendingLabels(const numSeqs_t capacity);

        /*
         * Build the data structure for the provided sequence of integers.
         */
        RankedAscendingLabels(std::vector<numSeqs_t>& labels);

        ~RankedAscendingLabels();

        RankedAscendingLabels(const RankedAscendingLabels& other); // copy constructor

        RankedAscendingLabels(RankedAscendingLabels&& other) noexcept; // move constructor

        RankedAscendingLabels& operator=(const RankedAscendingLabels& other); // copy assignment operator

        RankedAscendingLabels& operator=(RankedAscendingLabels&& other) noexcept; // move assignment operator

        RankedAscendingLabels* clone() const; // deep-copy clone method

        /*
         * Adds the next value to the data structure and returns it ranked value.
         * Does not check whether there is space for another value.
         */
        numSeqs_t add(const numSeqs_t lab);

        /*
         * Determine the original value of the ranked value.
         */
        numSeqs_t unrank(const numSeqs_t r) const;

        /*
         * Check whether the stored sequence contains the provided original value.
         */
        bool contains(const numSeqs_t lab) const;

        /*
         * Check whether the stored sequence contains the provided ranked value.
         */
        bool contains_rank(const numSeqs_t r) const;

        /*
         * Remove the provided original value from the sequence and return its rank.
         */
        numSeqs_t remove(const numSeqs_t lab);

        /*
         * Returns the size of the currently stored sequence of integers.
         */
        numSeqs_t size() const;

        /*
         * Determine the memory consumption (in bytes).
         */
        size_t size_in_bytes() const;

    private:
        long long* labels_; // the sequence of original integers
        numSeqs_t size_; // the current length of the stored sequence
        numSeqs_t capa_; // the maximum length of the stored sequence

    };

    /*
     * Collection of (inverted) indices for the segment filter.
     * Each inverted index is a binary relation of type R between segments (strings) and identifiers.
     *
     * To be used with k^2-trees (see K2Tree and the derived classes) for the binary relations.
     *
     * The collection of indices is built for an amplicon pool and segment filters using a given number of segments.
     * The provided operations are also tailored to the use with the implemented segment filters in AuxiliaryData.
     * Therefore, some operations expect a certain order w.r.t. the provided parameters.
     * For example, when adding entries the identifiers are expected to be monotonically increasing.
     */
    template<typename R>
    class SpaceSegmentRelations {

        typedef std::unordered_map<numSeqs_t, numSeqs_t> length_map_type;
        typedef std::unordered_map<std::unique_ptr<SubstringWrapper>, numSeqs_t, hashSubstringWrapper, equalSubstringWrapper> sip_map_type;

    public:
        typedef R rel_type;
        typedef typename rel_type::positions_t pairs_type;

        SpaceSegmentRelations(const AmpliconCollection& ac, numSeqs_t num_segs, const K2TreeParameters& parameters) {

            pool_size_ = ac.size();
            num_segs_ = num_segs;
            tree_parameters_ = parameters;

            std::set<lenSeqs_t> lengths;
            for (numSeqs_t i = 0; i < pool_size_; i++) {
                lengths.insert(ac.len(i));
            }

            numSeqs_t offset = 0;
            for (auto len : lengths) {
                access_len_[len] = offset;
                offset += num_segs_;
            }

            access_sip_ = new sip_map_type[access_len_.size() * num_segs_];
            relations_ = new rel_type[access_len_.size() * num_segs_];
            labels_ = new RankedAscendingLabels[lengths.size()];

        }

        virtual ~SpaceSegmentRelations() {

            delete[] labels_;
            delete[] relations_;
            delete[] access_sip_;

        }

        SpaceSegmentRelations(const SpaceSegmentRelations& other) : // copy constructor
                pool_size_(other.pool_size_), num_segs_(other.num_segs_), access_len_(other.access_len_), tree_parameters_(other.tree_parameters_) {

            access_sip_ = new sip_map_type[access_len_.size() * num_segs_];
            for (auto i = 0; i < access_len_.size() * num_segs_; i++) {
                for (auto& kv : other.access_sip_[i]) {
                    access_sip_[i][std::unique_ptr<SubstringWrapper>(kv.first->clone())] = kv.second;
                }
            }

            relations_ = new rel_type[pool_size_ * num_segs_];
            std::copy(other.relations_, other.relations_ + pool_size_ * num_segs_, relations_);

            labels_ = new RankedAscendingLabels[access_len_.size()];
            std::copy(other.labels_, other.labels_ + access_len_.size(), labels_);

        }

        SpaceSegmentRelations(SpaceSegmentRelations&& other) noexcept : // move constructor
                pool_size_(other.pool_size_), num_segs_(other.num_segs_), tree_parameters_(other.tree_parameters_), access_len_(std::move(other.access_len_)) {

            access_sip_ = other.access_sip_; other.access_sip_ = nullptr;
            relations_ = other.relations_; other.relations_ = nullptr;
            labels_ = other.labels_; other.labels_ = nullptr;

        }

        SpaceSegmentRelations& operator=(const SpaceSegmentRelations& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] relations_;
            delete[] access_sip_;

            // copy new resources
            pool_size_ = other.pool_size_;
            num_segs_ = other.num_segs_;
            tree_parameters_ = other.tree_parameters_;

            access_len_ = other.access_len_;

            access_sip_ = new sip_map_type[access_len_.size() * num_segs_];
            for (auto i = 0; i < access_len_.size() * num_segs_; i++) {
                for (auto& kv : other.access_sip_[i]) {
                    access_sip_[i][std::unique_ptr<SubstringWrapper>(kv.first->clone())] = kv.second;
                }
            }

            relations_ = new rel_type[pool_size_ * num_segs_];
            std::copy(other.relations_, other.relations_ + pool_size_ * num_segs_, relations_);

            labels_ = new RankedAscendingLabels[access_len_.size()];
            std::copy(other.labels_, other.labels_ + access_len_.size(), labels_);

            return *this;

        }

        SpaceSegmentRelations& operator=(SpaceSegmentRelations&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] relations_;
            delete[] access_sip_;

            // copy new resources
            pool_size_ = other.pool_size_;
            num_segs_ = other.num_segs_;
            tree_parameters_ = other.tree_parameters_;

            access_len_ = std::move(other.access_len_);
            access_sip_ = other.access_sip_; other.access_sip_ = nullptr;
            relations_ = other.relations_; other.relations_ = nullptr;
            labels_ = other.labels_; other.labels_ = nullptr;

            return *this;

        }

        // "count" the amplicons of the given length whose seg-th segment is equal to the given substring,
        // appends the ids of the found amplicons to the given vector
        // IMPORTANT: The existence of the length must be checked beforehand (or use add_counts_checked).
        virtual void add_counts(const numSeqs_t len, const numSeqs_t seg, const std::unique_ptr<SubstringWrapper>& sub_sip, std::vector<numSeqs_t>& cand_cnts) {

            auto& acc = access_sip_[access_len_[len] + seg];
            auto& labs = labels_[access_len_[len] / num_segs_];
            auto iter = acc.find(sub_sip);

            if (iter != acc.end()) {

                for (auto& s : relations_[access_len_[len] + seg].get_successors(iter->second)) {
                    cand_cnts.push_back(labs.unrank(s));
                }

            }

        }
        virtual void add_counts_checked(const numSeqs_t len, const numSeqs_t seg, const std::unique_ptr<SubstringWrapper>& sub_sip, std::vector<numSeqs_t>& cand_cnts) {

            auto& acc = access_sip_[access_len_[len] + seg];
            auto& labs = labels_[access_len_[len] / num_segs_];
            auto iter = acc.find(sub_sip);

            if (iter != acc.end()) {

                for (auto& s : relations_[access_len_[len] + seg].get_successors(iter->second)) {
                    cand_cnts.push_back(labs.unrank(s));
                }

            }

        }

        // "count" the amplicons of the given length whose seg-th segment is equal to the given substring
        // or the num_neigh next substrings generated by shifting to the right,
        // appends the ids of the found amplicons to the given vector
        // IMPORTANT: The existence of the length must be checked beforehand (or use add_counts_checked).
        virtual void add_neighboured_counts(const numSeqs_t len, const numSeqs_t seg, std::unique_ptr<SubstringWrapper>& sub_sip,
                                    numSeqs_t num_neigh, std::vector<numSeqs_t>& cand_cnts) {

            auto& acc = access_sip_[access_len_[len] + seg];
            auto& labs = labels_[access_len_[len] / num_segs_];

            for (auto n = 0; n <= num_neigh; n++) {

                auto iter = acc.find(sub_sip);

                if (iter != acc.end()) {

                    for (auto& s : relations_[access_len_[len] + seg].get_successors(iter->second)) {
                        cand_cnts.push_back(labs.unrank(s));
                    }

                }

                sub_sip->shift();

            }

        }
        virtual void add_neighboured_counts_checked(const numSeqs_t len, const numSeqs_t seg, std::unique_ptr<SubstringWrapper>& sip,
                                            numSeqs_t num_neigh, std::vector<numSeqs_t>& cand_cnts) {

            if (access_len_.find(len) != access_len_.end()) {

                auto& acc = access_sip_[access_len_[len] + seg];
                auto& labs = labels_[access_len_[len] / num_segs_];

                for (auto n = 0; n <= num_neigh; n++) {

                    auto iter = acc.find(sip);

                    if (iter != acc.end()) {

                        for (auto& s : relations_[access_len_[len] + seg].get_successors(iter->second)) {
                            cand_cnts.push_back(labs.unrank(s));
                        }

                    }

                    sip->shift();

                }

            }

        }

        // inserts the given segment in the relation for the given length and segment number
        numSeqs_t register_segment(const numSeqs_t len, const numSeqs_t seg, const std::unique_ptr<SubstringWrapper>& sub_sip) {

            auto& acc = access_sip_[access_len_[len] + seg];

            auto iter = acc.find(sub_sip);
            if (iter == acc.end()) {

                size_t s = acc.size();
                acc[std::unique_ptr<SubstringWrapper>(sub_sip->clone())] = s;
                return s;

            } else {
                return iter->second;
            }

        }

        void build_labels(const numSeqs_t len, std::vector<numSeqs_t>& labels) {
            labels_[access_len_[len] / num_segs_] = RankedAscendingLabels(labels);
        }

        void build_relation(const numSeqs_t len, const numSeqs_t seg, pairs_type& pairs) {
            relations_[access_len_[len] + seg] = rel_type(pairs, tree_parameters_);
        }

        // remove the specified amplicon from all relations,
        // expects that only previously inserted identifiers are removed exactly once
        // and that the identifiers were inserted in ascending order
        virtual void remove_amplicon(const lenSeqs_t len, const numSeqs_t ampl_id) {

            auto r = labels_[access_len_[len] / num_segs_].remove(ampl_id);

            for (numSeqs_t seg = 0; seg < num_segs_; seg++) {
                relations_[access_len_[len] + seg].set_null_column(r, tree_parameters_.deep);
            }

        }

        // check whether there are relations for the given length
        bool contains_length(const numSeqs_t len) const {
            return access_len_.find(len) != access_len_.end();
        }

        /*
         * Determine the memory consumption (in bytes).
         */
        virtual size_t size_in_bytes() const {

            size_t allocated_access_sip = 0;
            for (auto i = 0; i < access_len_.size() * num_segs_; i++) {
                allocated_access_sip += sip_map_single_size_in_bytes(access_sip_[i]);
            }
            size_t allocated_relations = 0;
            for (auto i = 0; i < access_len_.size() * num_segs_; i++) {
                allocated_relations += relations_[i].size_in_bytes();
            }
            size_t allocated_labels = 0;
            for (auto i = 0; i < access_len_.size(); i++) {
                allocated_labels += labels_[i].size_in_bytes();
            }

            return sizeof(SpaceSegmentRelations) // SpaceSegmentRelations itself
                + allocated_access_sip // access_sip_
                + allocated_relations // relations_
                + allocated_labels // labels_
                + umap_pod_size_in_bytes(access_len_) - sizeof(length_map_type); // access_len_ (do not count size of length_map_type object twice)

        }

        /*
         * Show memory profile.
         */
        virtual void show_memory() const {

            size_t allocated_access_sip = 0;
            for (auto i = 0; i < access_len_.size() * num_segs_; i++) {
                allocated_access_sip += sip_map_single_size_in_bytes(access_sip_[i]);
            }
            size_t allocated_relations = 0;
            for (auto i = 0; i < access_len_.size() * num_segs_; i++) {
                allocated_relations += relations_[i].size_in_bytes();
            }
            size_t allocated_labels = 0;
            for (auto i = 0; i < access_len_.size(); i++) {
                allocated_labels += labels_[i].size_in_bytes();
            }

            std::cout << "#  * (Static)SpaceSegmentRelations " << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * sizeof(SpaceSegmentRelations): " << sizeof(SpaceSegmentRelations) << " bytes" << std::endl;
            std::cout << "#  * sizeof(sip_map_type): " << sizeof(sip_map_type) << " bytes" << std::endl;
            std::cout << "#  * sizeof(sip_map_type*): " << sizeof(sip_map_type*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(rel_type): " << sizeof(rel_type) << " bytes" << std::endl;
            std::cout << "#  * sizeof(rel_type*): " << sizeof(rel_type*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(RankedAscendingLabels): " << sizeof(RankedAscendingLabels) << " bytes" << std::endl;
            std::cout << "#  * sizeof(RankedAscendingLabels*): " << sizeof(RankedAscendingLabels*) << " bytes" << std::endl;
            std::cout << "#  * sizeof(length_map_type): " << sizeof(length_map_type) << " bytes" << std::endl;
            std::cout << "#  * sizeof(numSeqs_t): " << sizeof(numSeqs_t) << " bytes" << std::endl;
            std::cout << "#  * sizeof(K2TreeParameters): " << sizeof(K2TreeParameters) << " bytes" << std::endl;
            std::cout << "#  * " << std::endl;
            std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
            std::cout << "#  * access_sip_: " << allocated_access_sip << " bytes" << std::endl;
            std::cout << "#  * relations_: " << allocated_relations << " bytes" << std::endl;
            std::cout << "#  * labels_: " << allocated_labels << " bytes" << std::endl;
            std::cout << "#  * access_len_: " << umap_pod_size_in_bytes(access_len_) << " bytes" << std::endl;

        }

    protected:
        // array of mappings from SequenceWrapper to (numeric) id,
        // contains a mapping for each observed sequence length and segment number,
        // mappings belonging to the same sequence length are stored consecutively
        sip_map_type* access_sip_;

        // array of binary relations,
        // contains a binary relation for each observed sequence length and segment number,
        // relations belonging to the same sequence length are stored consecutively
        rel_type* relations_;

        // array of mappings reducing the identifiers to consecutive ranges,
        // contains a mapping for each observed sequence length,
        // allow to record whether a label (and associated columns in the binary relations) have already been removed
        RankedAscendingLabels* labels_;

        // mapping from sequence length to the first entry in access_sip_ resp. relations_ belonging to that length
        length_map_type access_len_;

        numSeqs_t pool_size_; // number of amplicons in the pool for which this relation has been built
        numSeqs_t num_segs_; // number of segments per sequence

        K2TreeParameters tree_parameters_; // parameters of the k^2-trees for the binary relations

    };

    /*
     * Variant of above SpaceSegmentRelations.
     *
     * Works in the same way, except that it does not modify the binary relations (by removing columns)
     * but rather filters the returned identifiers using the RankedAscendingLabels before passing them onto the caller.
     */
    template<typename R>
    class StaticSpaceSegmentRelations : public SpaceSegmentRelations<R> {

        using SpaceSegmentRelations<R>::access_sip_;
        using SpaceSegmentRelations<R>::relations_;
        using SpaceSegmentRelations<R>::access_len_;
        using SpaceSegmentRelations<R>::pool_size_;
        using SpaceSegmentRelations<R>::num_segs_;
        using SpaceSegmentRelations<R>::labels_;
        using SpaceSegmentRelations<R>::tree_parameters_;

    public:
        StaticSpaceSegmentRelations(const AmpliconCollection& ac, numSeqs_t num_segs, const K2TreeParameters& parameters)
            : SpaceSegmentRelations<R>(ac, num_segs, parameters) {
            // nothing else to do
        }

        virtual ~StaticSpaceSegmentRelations() {
            // nothing to do
        }

        StaticSpaceSegmentRelations(const StaticSpaceSegmentRelations& other) : SpaceSegmentRelations<R>(other) { // copy constructor
            // nothing else to do
        }

        StaticSpaceSegmentRelations(StaticSpaceSegmentRelations&& other) noexcept : SpaceSegmentRelations<R>(other) { // move constructor
            // nothing else to do
        }

        StaticSpaceSegmentRelations& operator=(const StaticSpaceSegmentRelations& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            SpaceSegmentRelations<R>::operator=(other);

            return *this;

        }

        StaticSpaceSegmentRelations& operator=(StaticSpaceSegmentRelations&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            SpaceSegmentRelations<R>::operator=(other);

            return *this;

        }

        // "count" the amplicons of the given length whose seg-th segment is equal to the given substring,
        // appends the ids of the found amplicons to the given vector if they are not marked as removed
        // IMPORTANT: The existence of the length must be checked beforehand (or use add_counts_checked).
        void add_counts(const numSeqs_t len, const numSeqs_t seg, const std::unique_ptr<SubstringWrapper>& sub_sip,
                        std::vector<numSeqs_t>& cand_cnts) override {

            auto& acc = access_sip_[access_len_[len] + seg];
            auto iter = acc.find(sub_sip);

            if (iter != acc.end()) {

                auto& labs = labels_[access_len_[len] / num_segs_];
                for (auto& s : relations_[access_len_[len] + seg].get_successors(iter->second)) {
                    if (labs.contains_rank(s)) cand_cnts.push_back(labs.unrank(s));
                }

            }

        }
        void add_counts_checked(const numSeqs_t len, const numSeqs_t seg, const std::unique_ptr<SubstringWrapper>& sub_sip,
                                std::vector<numSeqs_t>& cand_cnts) override {

            auto& acc = access_sip_[access_len_[len] + seg];
            auto iter = acc.find(sub_sip);

            if (iter != acc.end()) {

                auto& labs = labels_[access_len_[len] / num_segs_];
                for (auto& s : relations_[access_len_[len] + seg].get_successors(iter->second)) {
                    if (labs.contains_rank(s)) cand_cnts.push_back(labs.unrank(s));
                }

            }

        }

        // "count" the amplicons of the given length whose seg-th segment is equal to the given substring
        // or the num_neigh next substrings generated by shifting to the right,
        // appends the ids of the found amplicons to the given vector if they are not marked as removed
        // IMPORTANT: The existence of the length must be checked beforehand (or use add_counts_checked).
        void add_neighboured_counts(const numSeqs_t len, const numSeqs_t seg, std::unique_ptr<SubstringWrapper>& sub_sip,
                                    numSeqs_t num_neigh, std::vector<numSeqs_t>& cand_cnts) override {

            auto& acc = access_sip_[access_len_[len] + seg];
            auto& labs = labels_[access_len_[len] / num_segs_];

            for (auto n = 0; n <= num_neigh; n++) {

                auto iter = acc.find(sub_sip);

                if (iter != acc.end()) {

                    for (auto& s : relations_[access_len_[len] + seg].get_successors(iter->second)) {
                        if (labs.contains_rank(s)) cand_cnts.push_back(labs.unrank(s));
                    }

                }

                sub_sip->shift();

            }

        }
        void add_neighboured_counts_checked(const numSeqs_t len, const numSeqs_t seg, std::unique_ptr<SubstringWrapper>& sip,
                                            numSeqs_t num_neigh, std::vector<numSeqs_t>& cand_cnts) override {

            if (access_len_.find(len) != access_len_.end()) {

                auto& acc = access_sip_[access_len_[len] + seg];
                auto& labs = labels_[access_len_[len] / num_segs_];

                for (auto n = 0; n <= num_neigh; n++) {

                    auto iter = acc.find(sip);

                    if (iter != acc.end()) {

                        for (auto& s : relations_[access_len_[len] + seg].get_successors(iter->second)) {
                            if (labs.contains_rank(s)) cand_cnts.push_back(labs.unrank(s));
                        }

                    }

                    sip->shift();

                }

            }

        }

        // mark the specified amplicon as removed from all relations,
        // expects that only previously inserted identifiers are removed exactly once
        // and that the identifiers were inserted in ascending order
        void remove_amplicon(const lenSeqs_t len, const numSeqs_t ampl_id) override {
            labels_[access_len_[len] / num_segs_].remove(ampl_id);
        }

    };


    /*
     * Support for clustering phases employing a (one-way) segment filter to efficiently search
     * for (not yet swarmed) partners of an amplicon.
     *
     * Similar to SegmentFilterAuxiliaryData but designed for the use with SpaceSegmentRelations
     * and StaticSpaceSegmentRelations.
     *
     * The implementation of the segment filter is based on:
     * Li et al. (2013), A partition-based method for string similarity joins with edit-distance constraints
     * Lin et al. (2014), Large-Scale Similarity Join with Edit-Distance Constraints
     *
     * Keeps recorded sequences in a set of strings.
     */
    template<typename S>
    class SpaceSegmentFilterAuxiliaryData : public AuxiliaryData {

        // maps sequence segments to amplicon 'ids' (represented by their indices within the AmpliconPool)
        typedef S SwarmingIndices;

    public:
        /*
         * Prepare the auxiliary data for the clustering phase (inverted indices on all amplicons).
         */
        SpaceSegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage, const numSeqs_t pool_id,
                                        lenSeqs_t threshold, lenSeqs_t num_extra_segments, const Configuration& config) :
                ac_(amplicon_storage.get_pool(pool_id)),
                indices_(ac_, threshold + num_extra_segments,
                         static_cast<const SpaceLevenshteinConfiguration&>(config).tree_parameters) {

            threshold_ = threshold;
            num_extra_segments_ = num_extra_segments;
            break_swarms_ = config.break_swarms;

            std::map<lenSeqs_t, std::vector<numSeqs_t>> length_groups;
            for (auto i = 0; i < ac_.size(); i++) {
                length_groups[ac_.len(i)].push_back(i);
            }

            prepare(length_groups);

            dist_fun_ = config.build_distance_function(amplicon_storage, threshold);

        }

        /*
         * Prepare the auxiliary data for the refinement phase (inverted indices on amplicon from light swarms).
         */
        SpaceSegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                                        const numSeqs_t pool_id, lenSeqs_t threshold, lenSeqs_t num_extra_segments,
                                        const Configuration& config) :
                ac_(amplicon_storage.get_pool(pool_id)),
                indices_(ac_, threshold + num_extra_segments,
                         static_cast<const SpaceLevenshteinConfiguration&>(config).tree_parameters) {

            threshold_ = threshold;
            num_extra_segments_ = num_extra_segments;
            break_swarms_ = false;

            prepare_light(swarm_storage.get_swarms(pool_id), config.boundary);

            dist_fun_ = config.build_distance_function(amplicon_storage, threshold);

        }

        virtual ~SpaceSegmentFilterAuxiliaryData() {

            delete dist_fun_;

            for (auto& kv : substrs_archive_) {
                delete[] kv.second; kv.second = nullptr;
            }

        }

        SpaceSegmentFilterAuxiliaryData(const SpaceSegmentFilterAuxiliaryData& other) : // copy constructor
                ac_(other.ac_), dist_fun_(other.dist_fun_->clone()),
                threshold_(other.threshold_), num_extra_segments_(other.num_extra_segments_),
                break_swarms_(other.break_swarms_), indices_(other.indices_) {

            for (auto& s : other.different_seqs_) different_seqs_.emplace(s->clone());

            for (auto& kv : other.substrs_archive_) {

                auto tmp = new Substrings[threshold_ + num_extra_segments_];
                std::copy(kv.second, kv.second + threshold_ + num_extra_segments_, tmp);
                substrs_archive_[kv.first] = tmp;

            }

        }

        SpaceSegmentFilterAuxiliaryData(SpaceSegmentFilterAuxiliaryData&& other) noexcept : // move constructor
                ac_(other.ac_), dist_fun_(other.dist_fun_), different_seqs_(std::move(other.different_seqs_)),
                threshold_(other.threshold_), num_extra_segments_(other.num_extra_segments_),
                break_swarms_(other.break_swarms_), indices_(std::move(other.indices_)),
                substrs_archive_(std::move(other.substrs_archive_)) {

            other.dist_fun_ = nullptr;

        }

        SpaceSegmentFilterAuxiliaryData& operator=(const SpaceSegmentFilterAuxiliaryData& other) = delete; // copy assignment operator

        SpaceSegmentFilterAuxiliaryData& operator=(SpaceSegmentFilterAuxiliaryData&& other) = delete; // move assignment operator

        SpaceSegmentFilterAuxiliaryData* clone() const override { // deep-copy clone method
            return new SpaceSegmentFilterAuxiliaryData(*this);
        }

        void tick_off_amplicon(numSeqs_t ampl_id) override {
            indices_.remove_amplicon(ac_.len(ampl_id), ampl_id);
        }

        bool record_amplicon(numSeqs_t ampl_id) override {
            return different_seqs_.emplace(ac_.seq(ampl_id)).second;
        }

        void clear_amplicon_records() override {
            different_seqs_.clear();
        }

        std::vector<Partner> find_partners(numSeqs_t ampl_id) override {

            std::vector<Partner> partners;
            std::vector<numSeqs_t> cand_cnts;

            auto ampl_seq = ac_.seq(ampl_id);
            numSeqs_t ampl_ab = ac_.ab(ampl_id);

            for (lenSeqs_t seq_len = ac_.len(ampl_id), partner_len = (seq_len > threshold_) * (seq_len - threshold_);
                 partner_len <= seq_len + threshold_; partner_len++) {

                cand_cnts.clear();

                add_candidate_counts(*ampl_seq, seq_len, partner_len, cand_cnts);
                verify_candidates(ampl_id, ampl_ab, cand_cnts, partners);

            }

            return partners;

        }

        size_t size_in_bytes() const override {
            return sizeof(SpaceSegmentFilterAuxiliaryData) // SpaceSegmentFilterAuxiliaryData itself
                + dist_fun_->size_in_bytes() // dist_fun_
                + different_seqs_size_in_bytes(different_seqs_) - sizeof(std::set<std::unique_ptr<SequenceWrapper>, lessSequenceWrapper>) // different_seqs_ (do not count size of set object twice)
                + indices_.size_in_bytes() - sizeof(SwarmingIndices) // indices_ (do not count size of SwarmingIndices object twice)
                + substrs_archive_size_in_bytes(substrs_archive_) - sizeof(std::map<std::pair<lenSeqs_t, lenSeqs_t>, Substrings*>); // substr_archive_ (do not count size of map object twice)
        }

        void show_memory(numSeqs_t pid) const {

            std::cout << "##################################################" << std::endl;
            std::cout << "# SpaceSegmentFilterAuxiliaryData " << std::endl;
            std::cout << "#  - pid: " << pid << std::endl;
            std::cout << "# " << std::endl;
            std::cout << "# sizeof(SpaceSegmentFilterAuxiliaryData): " << sizeof(SpaceSegmentFilterAuxiliaryData) << " bytes" << std::endl;
            std::cout << "# sizeof(Distance*): " << sizeof(Distance*) << " bytes" << std::endl;
            std::cout << "# sizeof(std::set<std::unique_ptr<SequenceWrapper>, lessSequenceWrapper>): " << sizeof(std::set<std::unique_ptr<SequenceWrapper>, lessSequenceWrapper>) << " bytes" << std::endl;
            std::cout << "# sizeof(lenSeqs_t): " << sizeof(lenSeqs_t) << " bytes" << std::endl;
            std::cout << "# sizeof(bool): " << sizeof(bool) << " bytes" << std::endl;
            std::cout << "# sizeof(SwarmingIndices): " << sizeof(SwarmingIndices) << " bytes" << std::endl;
            std::cout << "# sizeof(std::map<std::pair<lenSeqs_t, lenSeqs_t>, Substrings*>): " << sizeof(std::map<std::pair<lenSeqs_t, lenSeqs_t>, Substrings*>) << " bytes" << std::endl;
            std::cout << "# " << std::endl;
            std::cout << "# Total size: " << size_in_bytes() << " bytes" << std::endl;
            std::cout << "# dist_fun_: " << dist_fun_->size_in_bytes() << " bytes" << std::endl;
            std::cout << "# different_seqs_: " << different_seqs_size_in_bytes(different_seqs_) << " bytes" << std::endl;
            std::cout << "# indices_: " << indices_.size_in_bytes() << " bytes" << std::endl;
            indices_.show_memory();
            std::cout << "# substrs_archive_: " << substrs_archive_size_in_bytes(substrs_archive_) << " bytes" << std::endl;
            std::cout << "##################################################" << std::endl;

        }

    private:
        /*
         * Prepare inverted indices using all amplicons from the amplicon collection.
         */
        void prepare(std::map<lenSeqs_t, std::vector<numSeqs_t>>& length_groups) {

            // index all amplicons, grouped by sequence length
            for (auto& kv : length_groups) {

                lenSeqs_t seq_len = kv.first;
                auto& ids = kv.second;

                if (substrs_archive_.find(std::make_pair(seq_len, seq_len)) == substrs_archive_.end()) {

                    // substrings information
                    std::pair<lenSeqs_t, lenSeqs_t> len_pair(seq_len, 0);
                    for (lenSeqs_t partner_len = (seq_len > threshold_) * (seq_len - threshold_);
                         partner_len <= seq_len + threshold_; partner_len++) {

                        len_pair.second = partner_len;
                        auto arr = new Substrings[threshold_ + num_extra_segments_];
                        for (lenSeqs_t segment_index = 0; segment_index < threshold_ + num_extra_segments_; segment_index++) {
                            if (partner_len <= seq_len) {
                                arr[segment_index] = select_substrs(seq_len, partner_len, segment_index, threshold_, num_extra_segments_);
                            } else {
                                arr[segment_index] = select_substrs_backward(seq_len, partner_len, segment_index, threshold_, num_extra_segments_);
                            }
                        }
                        substrs_archive_[len_pair] = arr;

                    }

                }

                lenSeqs_t d = seq_len - (seq_len / (threshold_ + num_extra_segments_)) * (threshold_ + num_extra_segments_);
                lenSeqs_t l = (seq_len / (threshold_ + num_extra_segments_));

                indices_.build_labels(seq_len, ids);

                for (lenSeqs_t s = 0, start = 0, len = l + (s >= (threshold_ + num_extra_segments_ - d));
                     s < threshold_ + num_extra_segments_;
                     s++, start += len, len = l + (s >= (threshold_ + num_extra_segments_ - d))) {

                    typename SwarmingIndices::pairs_type pairs;
                    for (auto i = 0; i < ids.size(); i++) {
                        pairs.emplace_back(indices_.register_segment(seq_len, s, ac_.seq(ids[i])->substr(start, len)), i);
                    }
                    indices_.build_relation(seq_len, s, pairs);

                }

            }

        }

        /*
         * Prepare inverted indices using only amplicons from light swarms.
         */
        void prepare_light(const Swarms& swarms, const numSeqs_t boundary) {

            std::map<lenSeqs_t, std::vector<numSeqs_t>> light_length_groups;
            std::set<lenSeqs_t> observed_lengths;

            // prepare indexing of amplicons of light swarms
            for (numSeqs_t s = 0; s < swarms.size(); s++) {

                auto& swarm = swarms.get(s);

                if (swarm.mass() < boundary && !swarm.is_attached()) {

                    for (numSeqs_t i = 0; i < swarm.size(); i++) {
                        light_length_groups[ac_.len(swarm.member(i))].push_back(swarm.member(i));
                    }

                } else {

                    for (numSeqs_t i = 0; i < swarm.size(); i++) {

                        lenSeqs_t seq_len = ac_.len(swarm.member(i));

                        if (observed_lengths.insert(seq_len).second) {

                            // substrings information of amplicons of heavy swarms
                            std::pair<lenSeqs_t, lenSeqs_t> len_pair(seq_len, 0);
                            for (lenSeqs_t partner_len = (seq_len > threshold_) * (seq_len - threshold_);
                                 partner_len <= seq_len + threshold_; partner_len++) {

                                len_pair.second = partner_len;
                                auto arr = new Substrings[threshold_ + num_extra_segments_];
                                for (lenSeqs_t segment_index = 0; segment_index < threshold_ + num_extra_segments_; segment_index++) {
                                    if (partner_len <= seq_len) {
                                        arr[segment_index] = select_substrs(seq_len, partner_len, segment_index, threshold_, num_extra_segments_);
                                    } else {
                                        arr[segment_index] = select_substrs_backward(seq_len, partner_len, segment_index, threshold_, num_extra_segments_);
                                    }
                                }
                                substrs_archive_[len_pair] = arr;

                            }

                        }

                    }

                }

            }

            // sort "labels" (pool-internal integer ids) of amplicons from light swarms per length group
            for (auto& kv : light_length_groups) {
                std::sort(kv.second.begin(), kv.second.end());
            }

            // index all light amplicons, grouped by sequence length
            for (auto& kv : light_length_groups) {

                lenSeqs_t seq_len = kv.first;
                auto& ids = kv.second;

                lenSeqs_t d = seq_len - (seq_len / (threshold_ + num_extra_segments_)) * (threshold_ + num_extra_segments_);
                lenSeqs_t l = (seq_len / (threshold_ + num_extra_segments_));

                indices_.build_labels(seq_len, ids);

                for (lenSeqs_t s = 0, start = 0, len = l + (s >= (threshold_ + num_extra_segments_ - d));
                     s < threshold_ + num_extra_segments_;
                     s++, start += len, len = l + (s >= (threshold_ + num_extra_segments_ - d))) {

                    typename SwarmingIndices::pairs_type pairs;
                    for (auto i = 0; i < ids.size(); i++) {
                        pairs.emplace_back(indices_.register_segment(seq_len, s, ac_.seq(ids[i])->substr(start, len)), i);
                    }
                    indices_.build_relation(seq_len, s, pairs);

                }

            }

        }

        /*
         * Search for segment matches with indexed amplicons of a given length and
         * add the found matches to the candidate counts.
         */
        void add_candidate_counts(const SequenceWrapper& ampl_seq, const lenSeqs_t ampl_len, lenSeqs_t partner_len, std::vector<numSeqs_t>& cand_cnts) {

            if (indices_.contains_length(partner_len)) {

                auto arr = substrs_archive_[std::make_pair(ampl_len, partner_len)];

                for (lenSeqs_t i = 0; i < threshold_ + num_extra_segments_; i++) {

                    Substrings& subs = arr[i];
                    auto sub_sip = ampl_seq.substr(subs.first, subs.len);

                    indices_.add_neighboured_counts(partner_len, i, sub_sip, subs.last - subs.first, cand_cnts);

                }

            }

        }

        /*
         * Verify the candidates by searching for amplicons with sufficient segment matches
         * and then computing the distance using the distance function.
         */
        void verify_candidates(const numSeqs_t ampl_id, const numSeqs_t ampl_ab, std::vector<numSeqs_t>& cand_cnts, std::vector<Partner>& partners) {

            std::sort(cand_cnts.begin(), cand_cnts.end());

            // general pigeonhole principle: for being a candidate, at least k segments have to be matched
            lenSeqs_t cnt = 0;
            numSeqs_t prev_cand = (cand_cnts.empty()) ? 0 : cand_cnts.front();
            for (auto cand_id : cand_cnts) {

                if (prev_cand != cand_id) {

                    if ((cnt >= num_extra_segments_) && (!break_swarms_ || ampl_ab >= ac_.ab(prev_cand))) {

                        auto dist = dist_fun_->distance(ac_, ampl_id, prev_cand);

                        if (dist <= threshold_) {
                            partners.emplace_back(prev_cand, dist);
                        }

                    }

                    cnt = 1;
                    prev_cand = cand_id;


                } else {
                    cnt++;
                }

            }

            if ((cnt >= num_extra_segments_) && (!break_swarms_ || ampl_ab >= ac_.ab(prev_cand))) {

                auto dist = dist_fun_->distance(ac_, ampl_id, prev_cand);

                if (dist <= threshold_) {
                    partners.emplace_back(prev_cand, dist);
                }

            }

        }


        const AmpliconCollection& ac_; // amplicon collection which is supported
        Distance* dist_fun_; // distance function used to find partners
        std::set<std::unique_ptr<SequenceWrapper>, lessSequenceWrapper> different_seqs_; // recorded different sequences

        lenSeqs_t threshold_; // threshold applied to partner search
        lenSeqs_t num_extra_segments_; // number of extra segments for segment filter
        bool break_swarms_; // flag indicating whether swarm breaking should be applied

        SwarmingIndices indices_; // inverted indices of the segment filter
        std::map<std::pair<lenSeqs_t, lenSeqs_t>, Substrings*> substrs_archive_; // substrings archive for more efficient segment filter

    };

}

#endif //GEFAST_FLEXIBLERELATIONS_HPP
