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

#ifndef GEFAST_AUXILIARYDATA_HPP
#define GEFAST_AUXILIARYDATA_HPP

#include "modes/AlignmentFreeMode.hpp"
#include "Base.hpp"
#include "Distances.hpp"
#include "Relations.hpp"

namespace GeFaST {

    /*
     * Simple support for realising the clustering phase.
     *
     * Keeps recorded sequences in a set of strings and maintains
     * a set of amplicons not yet included in a swarm.
     */
    class SimpleAuxiliaryData : public AuxiliaryData {

    public:
        SimpleAuxiliaryData(const AmpliconStorage& amplicon_storage, const numSeqs_t pool_id, const dist_t threshold,
                const Configuration& config);

        SimpleAuxiliaryData(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                const numSeqs_t pool_id, const dist_t threshold, const Configuration& config);

        virtual ~SimpleAuxiliaryData();

        SimpleAuxiliaryData(const SimpleAuxiliaryData& other); // copy constructor

        SimpleAuxiliaryData(SimpleAuxiliaryData&& other) noexcept; // move constructor

        SimpleAuxiliaryData& operator=(const SimpleAuxiliaryData& other) = delete; // copy assignment operator

        SimpleAuxiliaryData& operator=(SimpleAuxiliaryData&& other) = delete; // move assignment operator

        SimpleAuxiliaryData* clone() const override; // deep-copy clone method

        void tick_off_amplicon(numSeqs_t ampl_id) override;

        bool record_amplicon(numSeqs_t ampl_id) override;

        void clear_amplicon_records() override;

        std::vector<Partner> find_partners(numSeqs_t ampl_id) override;

        size_t size_in_bytes() const override;

        void show_memory(numSeqs_t pid) const override;

    protected:
        const AmpliconCollection& ac_; // amplicon collection which is supported
        Distance* dist_fun_; // distance function used to find partners
        std::set<std::unique_ptr<SequenceWrapper>, lessSequenceWrapper> different_seqs_; // recorded different sequences

        dist_t threshold_; // threshold applied during search for partners / similar amplicons
        bool break_swarms_; // flag indicating whether swarm breaking should be applied

        std::set<numSeqs_t> unswarmed_; // pool-internal integer ids of amplicons not yet included in a swarm

    };



    /*
     * Set of substrings chosen for comparison with a segment.
     */
    struct Substrings {

        lenSeqs_t first = 0; // start position of first substring to be checked
        lenSeqs_t last = 0; // start position of last substring to be checked
        lenSeqs_t len = 0; // common length of all substrings to be checked

        Substrings() = default;

        Substrings(lenSeqs_t fp, lenSeqs_t lp, lenSeqs_t l);

        /*
         * Determine the memory consumption (in bytes).
         */
        size_t size_in_bytes() const {
            return sizeof(Substrings);
        }

    };

    /*
     * Select 'substrings' for segment filter using the  MMASS (multimatch-aware substring selection) method
     * proposed in "A partition-based method for string similarity joins with edit-distance constraints"
     * (Li et al., 2013, https://doi.org/10.1145/2487259.2487261).
     * Reduces the number of selected substrings by considering the length, relative shift and the implied
     * existence of further matching substrings.
     */
    Substrings select_substrs(const lenSeqs_t self_len, const lenSeqs_t partner_len, const lenSeqs_t seg_index,
            const lenSeqs_t t, const lenSeqs_t k);
    Substrings select_substrs_backward(const lenSeqs_t self_len, const lenSeqs_t partner_len, const lenSeqs_t seg_index,
            const lenSeqs_t t, const lenSeqs_t k);

    /*
     * Approximate the memory consumption (in bytes) of an std::map<std::pair<lenSeqs_t, lenSeqs_t>, Substrings*> instance
     * as used by auxiliary data structures.
     * Assumes the _Rb_tree implementation of GCC 4.9.2.
     *
     * Inspired by https://stackoverflow.com/questions/720507/how-can-i-estimate-memory-usage-of-stdmap/720520.
     */
    size_t substrs_archive_size_in_bytes(const std::map<std::pair<lenSeqs_t, lenSeqs_t>, Substrings*>& map);


    /*
     * Set of segments to be indexed as (first position, length of segment)
     */
    typedef std::vector<std::pair<lenSeqs_t, lenSeqs_t>> Segments;

    /*
     * Select 'segments' for indexing step (are returned via the reference parameter).
     */
    void select_segments(Segments& segments, const lenSeqs_t seq_len, const lenSeqs_t t, const lenSeqs_t k);


    /*
     * Support for clustering phases employing a (one-way) segment filter to efficiently search
     * for (not yet swarmed) partners of an amplicon.
     *
     * The implementation of the segment filter is based on:
     * Li et al. (2013), A partition-based method for string similarity joins with edit-distance constraints
     * Lin et al. (2014), Large-Scale Similarity Join with Edit-Distance Constraints
     *
     * Keeps recorded sequences in a set of strings.
     */
    class SegmentFilterAuxiliaryData : public AuxiliaryData {

        // maps sequence segments to amplicon 'ids' (represented by their indices within the AmpliconPool)
        typedef SegmentRelations SwarmingIndices;

    public:
        /*
         * Prepare the auxiliary data for the clustering phase (inverted indices on all amplicons).
         */
        SegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage, const numSeqs_t pool_id, lenSeqs_t threshold,
                lenSeqs_t num_extra_segments, const Configuration& config);

        /*
         * Prepare the auxiliary data for the refinement phase (inverted indices on amplicon from light swarms).
         */
        SegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                const numSeqs_t pool_id, lenSeqs_t threshold, lenSeqs_t num_extra_segments, const Configuration& config);

        virtual ~SegmentFilterAuxiliaryData();

        SegmentFilterAuxiliaryData(const SegmentFilterAuxiliaryData& other); // copy constructor

        SegmentFilterAuxiliaryData(SegmentFilterAuxiliaryData&& other) noexcept; // move constructor

        SegmentFilterAuxiliaryData& operator=(const SegmentFilterAuxiliaryData& other) = delete; // copy assignment operator

        SegmentFilterAuxiliaryData& operator=(SegmentFilterAuxiliaryData&& other) = delete; // move assignment operator

        SegmentFilterAuxiliaryData* clone() const override; // deep-copy clone method

        void tick_off_amplicon(numSeqs_t ampl_id) override;

        bool record_amplicon(numSeqs_t ampl_id) override;

        void clear_amplicon_records() override;

        std::vector<Partner> find_partners(numSeqs_t ampl_id) override;

        size_t size_in_bytes() const override;

        void show_memory(numSeqs_t pid) const override;

    private:
        /*
         * Prepare inverted indices using all amplicons from the amplicon collection.
         */
        void prepare();

        /*
         * Prepare inverted indices using only amplicons from light swarms.
         */
        void prepare_light(const Swarms& swarms, const numSeqs_t boundary);

        /*
         * Search for segment matches with indexed amplicons of a given length and
         * add the found matches to the candidate counts.
         */
        void add_candidate_counts(const SequenceWrapper& ampl_seq, const lenSeqs_t ampl_len, lenSeqs_t partner_len, std::vector<numSeqs_t>& cand_cnts);

        /*
         * Verify the candidates by searching for amplicons with sufficient segment matches
         * and then computing the distance using the distance function.
         */
        void verify_candidates(const numSeqs_t ampl_id, const numSeqs_t ampl_ab, std::vector<numSeqs_t>& cand_cnts, std::vector<Partner>& partners);


        const AmpliconCollection& ac_; // amplicon collection which is supported
        Distance* dist_fun_; // distance function used to find partners
        std::set<std::unique_ptr<SequenceWrapper>, lessSequenceWrapper> different_seqs_; // recorded different sequences

        lenSeqs_t threshold_; // threshold applied to partner search
        lenSeqs_t num_extra_segments_; // number of extra segments for segment filter
        bool break_swarms_; // flag indicating whether swarm breaking should be applied

        SwarmingIndices indices_; // inverted indices of the segment filter
        std::map<std::pair<lenSeqs_t, lenSeqs_t>, Substrings*> substrs_archive_; // substrings archive for more efficient segment filter

    };


    /*
     * Support for clustering phases employing a (two-way) segment filter to efficiently search
     * for (not yet swarmed) partners of an amplicon.
     *
     * The implementation of the segment filter is based on:
     * Li et al. (2013), A partition-based method for string similarity joins with edit-distance constraints
     * Lin et al. (2014), Large-Scale Similarity Join with Edit-Distance Constraints
     * Huang et al. (2015), A Partition-Based Bi-directional Filtering Method for String Similarity JOINs
     *
     * Keeps recorded sequences in a set of strings.
     */
    class TwoWaySegmentFilterAuxiliaryData : public AuxiliaryData {

        // maps sequence segments to amplicon 'ids' (represented by their indices within the AmpliconPool)
        typedef SegmentRelations SwarmingIndices;

    public:
        /*
         * Prepare the auxiliary data for the clustering phase (inverted indices on all amplicons).
         */
        TwoWaySegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage, const numSeqs_t pool_id,
                lenSeqs_t threshold, lenSeqs_t num_extra_segments, const Configuration& config);

        /*
         * Prepare the auxiliary data for the refinement phase (inverted indices on amplicon from light swarms).
         */
        TwoWaySegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                const numSeqs_t pool_id, lenSeqs_t threshold, lenSeqs_t num_extra_segments, const Configuration& config);

        virtual ~TwoWaySegmentFilterAuxiliaryData();

        TwoWaySegmentFilterAuxiliaryData(const TwoWaySegmentFilterAuxiliaryData& other); // copy constructor

        TwoWaySegmentFilterAuxiliaryData(TwoWaySegmentFilterAuxiliaryData&& other) noexcept; // move constructor

        TwoWaySegmentFilterAuxiliaryData& operator=(const TwoWaySegmentFilterAuxiliaryData& other) = delete; // copy assignment operator

        TwoWaySegmentFilterAuxiliaryData& operator=(TwoWaySegmentFilterAuxiliaryData&& other) = delete; // move assignment operator

        TwoWaySegmentFilterAuxiliaryData* clone() const override; // deep-copy clone method

        void tick_off_amplicon(numSeqs_t ampl_id) override;

        bool record_amplicon(numSeqs_t ampl_id) override;

        void clear_amplicon_records() override;

        std::vector<Partner> find_partners(numSeqs_t ampl_id) override;

        size_t size_in_bytes() const override;

        void show_memory(numSeqs_t pid) const override;


    private:
        /*
         * Prepare inverted indices using all amplicons from the amplicon collection.
         */
        void prepare();

        /*
         * Prepare inverted indices using only amplicons from light swarms.
         */
        void prepare_light(const Swarms& swarms, const numSeqs_t boundary);

        /*
         * Search for segment matches with indexed amplicons of a given length and
         * add the found matches to the candidate counts.
         */
        void add_candidate_counts(const SequenceWrapper& ampl_seq, const lenSeqs_t ampl_len, lenSeqs_t partner_len, std::vector<numSeqs_t>& cand_cnts);

        /*
         * Verify the candidates by searching for amplicons with sufficient segment matches,
         * applying the segment filter in the other direction as well,
         * and then computing the distance using the distance function on the remaining candidates.
         */
        void verify_candidates(const numSeqs_t ampl_id, const lenSeqs_t ampl_len, const numSeqs_t ampl_ab,
                               const lenSeqs_t partner_len, std::vector<numSeqs_t>& cand_cnts, std::vector<Partner>& partners);


        const AmpliconCollection& ac_; // amplicon collection which is supported
        Distance* dist_fun_; // distance function used to find partners
        std::set<std::unique_ptr<SequenceWrapper>, lessSequenceWrapper> different_seqs_; // recorded different sequences

        lenSeqs_t threshold_; // threshold applied to partner search
        lenSeqs_t num_extra_segments_; // number of extra segments for segment filter
        bool break_swarms_; // flag indicating whether swarm breaking should be applied

        SwarmingIndices indices_; // inverted indices of the segment filter
        std::map<std::pair<lenSeqs_t, lenSeqs_t>, Substrings*> substrs_archive_; // substrings archive for more efficient segment filter

        Segments segments_; // segment information of currently searched amplicon
        std::vector<std::string> segment_strs_; // segment strings of currently searched amplicon

    };


    /*
     * Support for clustering phases employing a (one-way) segment filter to efficiently search
     * for (not yet swarmed) partners of an amplicon.
     *
     * The implementation of the segment filter is based on:
     * Li et al. (2013), A partition-based method for string similarity joins with edit-distance constraints
     * Lin et al. (2014), Large-Scale Similarity Join with Edit-Distance Constraints
     *
     * The number of segments depends on the number of edit operations in an alignment,
     * but - for this kind of auxiliary data - we assume that the distance between amplicons
     * is the score of the alignment (not the number of edit operations).
     * The scoring functions minimise the score (score of 0 for match, positive scores otherwise).
     * In order to determine the number of segments needed for the filter,
     * we bound the number of possible edit operations by how often the score of the 'cheapest'
     * edit operations fits into the given threshold.
     *
     * Keeps recorded sequences in a set of strings.
     */
    class ScoreSegmentFilterAuxiliaryData : public AuxiliaryData {

        // maps sequence segments to amplicon 'ids' (represented by their indices within the AmpliconPool)
        typedef SegmentRelations SwarmingIndices;

    public:
        /*
         * Prepare the auxiliary data for the clustering phase (inverted indices on all amplicons).
         */
        ScoreSegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage, const numSeqs_t pool_id,
                lenSeqs_t threshold, lenSeqs_t num_threshold_segments, lenSeqs_t num_extra_segments, const Configuration& config);

        /*
         * Prepare the auxiliary data for the refinement phase (inverted indices on amplicon from light swarms).
         */
        ScoreSegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                const numSeqs_t pool_id, lenSeqs_t threshold, lenSeqs_t num_threshold_segments,
                lenSeqs_t num_extra_segments, const Configuration& config);

        virtual ~ScoreSegmentFilterAuxiliaryData();

        ScoreSegmentFilterAuxiliaryData(const ScoreSegmentFilterAuxiliaryData& other); // copy constructor

        ScoreSegmentFilterAuxiliaryData(ScoreSegmentFilterAuxiliaryData&& other) noexcept; // move constructor

        ScoreSegmentFilterAuxiliaryData& operator=(const ScoreSegmentFilterAuxiliaryData& other) = delete; // copy assignment operator

        ScoreSegmentFilterAuxiliaryData& operator=(ScoreSegmentFilterAuxiliaryData&& other) = delete; // move assignment operator

        ScoreSegmentFilterAuxiliaryData* clone() const override; // deep-copy clone method


        void tick_off_amplicon(numSeqs_t ampl_id) override;

        bool record_amplicon(numSeqs_t ampl_id) override;

        void clear_amplicon_records() override;

        std::vector<Partner> find_partners(numSeqs_t ampl_id) override;

        size_t size_in_bytes() const override;

        void show_memory(numSeqs_t pid) const override;

    private:
        /*
         * Prepare inverted indices using all amplicons from the amplicon collection.
         */
        void prepare();

        /*
         * Prepare inverted indices using only amplicons from light swarms.
         */
        void prepare_light(const Swarms& swarms, const numSeqs_t boundary);

        /*
         * Search for segment matches with indexed amplicons of a given length and
         * add the found matches to the candidate counts.
         */
        void add_candidate_counts(const SequenceWrapper& ampl_seq, const lenSeqs_t ampl_len, lenSeqs_t partner_len, std::vector<numSeqs_t>& cand_cnts);

        /*
         * Verify the candidates by searching for amplicons with sufficient segment matches
         * and then computing the distance using the distance function.
         */
        void verify_candidates(const numSeqs_t ampl_id, const numSeqs_t ampl_ab, std::vector<numSeqs_t>& cand_cnts, std::vector<Partner>& partners);


        const AmpliconCollection& ac_; // amplicon collection which is supported
        Distance* dist_fun_; // distance function used to find partners
        std::set<std::unique_ptr<SequenceWrapper>, lessSequenceWrapper> different_seqs_; // recorded different sequences

        lenSeqs_t threshold_; // threshold applied to partner search
        lenSeqs_t num_threshold_segments_; // number of segments based on threshold and minimum penalty for edit operation for segment filter
        lenSeqs_t num_extra_segments_; // number of extra segments for segment filter
        bool break_swarms_; // flag indicating whether swarm breaking should be applied

        SwarmingIndices indices_; // inverted indices of the segment filter
        std::map<std::pair<lenSeqs_t, lenSeqs_t>, Substrings*> substrs_archive_; // substrings archive for more efficient segment filter

    };


    /*
     * Simple support for realising the clustering phase when working with features.
     *
     * Keeps recorded sequences in a set of strings and maintains
     * a set of amplicons not yet included in a swarm.
     *
     */
    class FeatureAuxiliaryData : public AuxiliaryData {

    public:
        FeatureAuxiliaryData(const FeatureAmpliconStorage& amplicon_storage, const numSeqs_t pool_id, const dist_t threshold,
                             const AlignmentFreeConfiguration& config);

        FeatureAuxiliaryData(const FeatureAmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                             const numSeqs_t pool_id, const dist_t threshold, const AlignmentFreeConfiguration& config);

        virtual ~FeatureAuxiliaryData();

        FeatureAuxiliaryData(const FeatureAuxiliaryData& other); // copy constructor

        FeatureAuxiliaryData(FeatureAuxiliaryData&& other) noexcept; // move constructor

        FeatureAuxiliaryData& operator=(const FeatureAuxiliaryData& other) = delete; // copy assignment operator

        FeatureAuxiliaryData& operator=(FeatureAuxiliaryData&& other) = delete; // move assignment operator

        FeatureAuxiliaryData* clone() const override; // deep-copy clone method

        void tick_off_amplicon(numSeqs_t ampl_id) override;

        bool record_amplicon(numSeqs_t ampl_id) override;

        void clear_amplicon_records() override;

        std::vector<Partner> find_partners(numSeqs_t ampl_id) override;

        size_t size_in_bytes() const override;

        void show_memory(numSeqs_t pid) const override;

    protected:
        const FeatureAmpliconCollection& ac_; // amplicon collection which is supported
        Distance* dist_fun_; // distance function used to find partners
        std::set<std::unique_ptr<SequenceWrapper>, lessSequenceWrapper> different_seqs_; // recorded different sequences

        dist_t threshold_; // threshold applied during search for partners / similar amplicons
        bool break_swarms_; // flag indicating whether swarm breaking should be applied

        std::set<numSeqs_t> unswarmed_; // pool-internal integer ids of amplicons not yet included in a swarm

    };


    /*
     * Support for clustering phases employing a space-partitioning structure to efficiently search
     * for (not yet swarmed) partners of an amplicon using the feature representations of the amplicons.
     *
     * Keeps recorded sequences in a set of strings.
     */
    template<typename T>
    class SpacePartitioningAuxiliaryData : public AuxiliaryData {

    public:
        /*
         * Prepare the auxiliary data for the clustering phase (space-partitioning structure on all amplicons).
         */
        SpacePartitioningAuxiliaryData(const FeatureAmpliconStorage& amplicon_storage, const numSeqs_t pool_id,
                                       const dist_t threshold, const AlignmentFreeConfiguration& config)
                : ac_(amplicon_storage.get_pool(pool_id)), swarmed_(ac_.size(), false) {

            threshold_ = threshold;
            break_swarms_ = config.break_swarms;
            dist_fun_ = config.build_distance_function(amplicon_storage, threshold);
            sps_ = new T(amplicon_storage.get_pool(pool_id));

        }

        /*
         * Prepare the auxiliary data for the refinement phase (space-partitioning structure on amplicons from light swarms).
         */
        SpacePartitioningAuxiliaryData(const FeatureAmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                                       const numSeqs_t pool_id, const dist_t threshold, const AlignmentFreeConfiguration& config) :
                ac_(amplicon_storage.get_pool(pool_id)), swarmed_(ac_.size(), false) {

            threshold_ = threshold;
            break_swarms_ = false;
            dist_fun_ = config.build_distance_function(amplicon_storage, threshold);

            std::vector<numSeqs_t> light_amplicons;
            auto& swarms = swarm_storage.get_swarms(pool_id);
            for (numSeqs_t s = 0; s < swarms.size(); s++) {

                auto& swarm = swarms.get(s);
                if (swarm.mass() < config.boundary && !swarm.is_attached()) { // check whether attached in case refinement is iterated

                    for (numSeqs_t i = 0; i < swarm.size(); i++) {
                        light_amplicons.push_back(swarm.member(i));
                    }

                }
            }

            sps_ = new T(amplicon_storage.get_pool(pool_id), light_amplicons);

        }

        ~SpacePartitioningAuxiliaryData() {

            delete sps_;
            delete dist_fun_;

        }

        SpacePartitioningAuxiliaryData(const SpacePartitioningAuxiliaryData& other) : ac_(other.ac_) { // copy constructor

            dist_fun_ = other.dist_fun_->clone();
            for (auto& s : other.different_seqs_) different_seqs_.emplace(s->clone());
            threshold_ = other.threshold_;
            break_swarms_ = other.break_swarms_;

            swarmed_ = other.swarmed_;
            sps_ = other.sps_->clone();


        }

        SpacePartitioningAuxiliaryData(SpacePartitioningAuxiliaryData&& other) noexcept : ac_(other.ac_), different_seqs_(std::move(other.different_seqs_)) { // move constructor

            dist_fun_ = other.dist_fun_; other.dist_fun_ = nullptr;
            threshold_ = other.threshold_;
            break_swarms_ = other.break_swarms_;

            swarmed_ = other.swarmed_;
            sps_ = other.sps_; other.sps_ = nullptr;

        }

        SpacePartitioningAuxiliaryData& operator=(const SpacePartitioningAuxiliaryData& other) = delete; // copy assignment operator

        SpacePartitioningAuxiliaryData& operator=(SpacePartitioningAuxiliaryData&& other) = delete; // move assignment operator

        SpacePartitioningAuxiliaryData* clone() const override { // deep-copy clone method
            return new SpacePartitioningAuxiliaryData(*this);
        }

        void tick_off_amplicon(numSeqs_t ampl_id) override {
            swarmed_[ampl_id] = true;
        }

        bool record_amplicon(numSeqs_t ampl_id) override {
            return different_seqs_.emplace(ac_.seq(ampl_id)).second;
        }

        void clear_amplicon_records() override {
            different_seqs_.clear();
        }

        std::vector<Partner> find_partners(numSeqs_t ampl_id) override {

            std::vector<Partner> partners = sps_->range_partner(ampl_id, threshold_);

            // exclude amplicons that are already part of a swarm
            auto end_keep = std::remove_if(partners.begin(), partners.end(), [this](Partner& p) {return swarmed_[p.id];});

            if (break_swarms_) { // also exclude amplicons with a higher abundance (if requested)

                numSeqs_t max_abund = ac_.ab(ampl_id);
                end_keep = std::remove_if(partners.begin(), end_keep, [this, max_abund](Partner& p) {return ac_.ab(p.id) > max_abund;});

            }

            partners.erase(end_keep, partners.end());

            return partners;

        }

        size_t size_in_bytes() const override {

            std::cerr << "ERROR: SpacePartitioningAuxiliaryData is currently not part of the space-efficiency analysis." << std::endl;
            return 0;

        }

        void show_memory(numSeqs_t pid) const override {
            std::cerr << "ERROR: SpacePartitioningAuxiliaryData is currently not part of the space-efficiency analysis." << std::endl;
        }

    protected:
        const FeatureAmpliconCollection& ac_; // amplicon collection which is supported
        Distance* dist_fun_; // distance function used to find partners
        std::set<std::unique_ptr<SequenceWrapper>, lessSequenceWrapper> different_seqs_; // recorded different sequences

        dist_t threshold_; // threshold applied to partner search
        bool break_swarms_; // flag indicating whether swarm breaking should be applied

        std::vector<bool> swarmed_; // flags indicating which amplicons are already included in a swarm
        T* sps_; // spatial partitioning structure used to efficiently query for partners

    };

}

#endif //GEFAST_AUXILIARYDATA_HPP
