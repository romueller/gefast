/*
 * GeFaST
 *
 * Copyright (C) 2016 - 2020 Robert Mueller
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

    protected:
        const AmpliconCollection& ac_; // amplicon collection which is supported
        Distance* dist_fun_; // distance function used to find partners
        std::set<std::string> different_seqs_; // recorded different sequences

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
     * Set of segments to be indexed as (first position, length of segment)
     */
    typedef std::vector<std::pair<lenSeqs_t, lenSeqs_t>> Segments;

    /*
     * Select 'segments' for indexing step (are returned via the reference parameter).
     */
    void select_segments(Segments& segments, const lenSeqs_t seq_len, const lenSeqs_t t, const lenSeqs_t k);


    /*
     * Pair of pointers (first, second) describing the (sub)string corresponding to [first, last).
     */
    typedef std::pair<const char*, const char*> StringIteratorPair;

    // hash function for StringIteratorPair
    struct hashStringIteratorPair {
        size_t operator()(const StringIteratorPair& p) const;
    };

    // comparison function for string equality
    struct equalStringIteratorPair {
        bool operator()(const StringIteratorPair& lhs, const StringIteratorPair& rhs) const;
    };

    // comparison function for lexicographical order
    struct lessStringIteratorPair {
        bool operator()(const StringIteratorPair& a, const StringIteratorPair& b) const;
    };


    /*
     * Collection of (inverted) indices for the segment filter.
     *
     * The inverted indices are arranged in a grid where
     * the columns correspond to segments and
     * the rows correspond to sequence lengths.
     *
     * During the execution of the segment filter,
     * 'older' rows can be removed once they correspond to
     * sequences too short to be able to provide candidates for
     * the current (and future) sequences.
     *
     */
    template<typename T>
    class RollingIndices {

    public:
        typedef std::vector<T> Row;


        RollingIndices(lenSeqs_t t, lenSeqs_t w, bool f, bool s = true) {

            threshold_ = t;
            width_ = w;
            forward_ = f;
            shrink_ = s;
            min_length_ = 0;
            max_length_ = 0;

            empty_ = T();
            empty_row_ = Row(0);

        }


        /*
         * Return the indices for the specified length.
         */
        Row& get_indices_row(const lenSeqs_t len) {

            auto iter = indices_.find(len);

            return (iter != indices_.end()) ? iter->second : empty_row_;

        }

        /*
         * Return the index corresponding to the specified length and segment
         */
        T& get_index(const lenSeqs_t len, const lenSeqs_t i) {

            if (i >= width_) return empty_;

            auto iter = indices_.find(len);

            return (iter != indices_.end()) ? (iter->second)[i] : empty_;

        }

        /*
         * Add new row (and remove then outdated rows).
         */
        void roll(const lenSeqs_t len) {

            if (indices_.find(len) == indices_.end()) {

                indices_[len] = Row(width_);
                min_length_ = std::min(min_length_, len);
                max_length_ = std::max(max_length_, len);
                if (shrink_) shrink(len);

            }

        }

        /*
         * Remove outdated rows.
         */
        void shrink(const lenSeqs_t cur) {

            if (forward_) {

                for (lenSeqs_t len = min_length_; len < (cur - threshold_); len++) {
                    indices_.erase(len);
                }
                for (lenSeqs_t len = cur - threshold_; len <= cur; len++) {

                    if (indices_.find(len) != indices_.end()) {

                        min_length_ = len;
                        break;

                    }

                }

            } else {

                for (lenSeqs_t len = cur + threshold_ + 1; len <= max_length_; len++) {
                    indices_.erase(len);
                }
                for (lenSeqs_t len = cur + threshold_; len >= cur; len--) {

                    if (indices_.find(len) != indices_.end()) {

                        max_length_ = len;
                        break;

                    }

                }

            }

        }

        /*
         * Check whether there are indices for the given length.
         */
        bool contains(const lenSeqs_t len) {
            return (indices_.find(len) != indices_.end());
        }

        /*
         * Return minimum length for which there are indices.
         */
        lenSeqs_t min_length() const {
            return min_length_;
        }

        /*
         * Return maximum length for which there are indices.
         */
        lenSeqs_t max_length() const {
            return max_length_;
        }

    private:
        lenSeqs_t threshold_; // limits number of rows when applying shrink()
        lenSeqs_t width_; // number of columns / segments per row
        lenSeqs_t min_length_; // minimum length for which there are indices
        lenSeqs_t max_length_; // maximum length for which there are indices

        std::unordered_map<lenSeqs_t, Row> indices_; // indices grid

        T empty_; // empty (dummy) index returned for out-of-bounds queries
        Row empty_row_; // empty (dummy) row returned for out-of-bounds queries

        bool forward_; // indicates whether rolling forwards (lengths increase) or backwards (lengths decrease)
        bool shrink_; // flag indicating whether roll() automatically shrinks the index

    };


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
        typedef LabelLinkBinaryRelation<StringIteratorPair, numSeqs_t, hashStringIteratorPair, equalStringIteratorPair> SwarmingInvertedIndex;
        typedef RollingIndices<SwarmingInvertedIndex> SwarmingIndices;

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
        void add_candidate_counts(const numSeqs_t ampl_id, lenSeqs_t partner_len, std::vector<numSeqs_t>& cand_cnts);

        /*
         * Verify the candidates by searching for amplicons with sufficient segment matches
         * and then computing the distance using the distance function.
         */
        void verify_candidates(const numSeqs_t ampl_id, std::vector<numSeqs_t>& cand_cnts, std::vector<Partner>& partners);


        const AmpliconCollection& ac_; // amplicon collection which is supported
        Distance* dist_fun_; // distance function used to find partners
        std::set<std::string> different_seqs_; // recorded different sequences

        lenSeqs_t threshold_; // threshold applied to partner search
        lenSeqs_t num_extra_segments_; // number of extra segments for segment filter
        bool break_swarms_; // flag indicating whether swarm breaking should be applied

        SwarmingIndices indices_; // inverted indices of the segment filter
        std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>> substrs_archive_; // substrings archive for more efficient segment filter

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
        typedef LabelLinkBinaryRelation<StringIteratorPair, numSeqs_t, hashStringIteratorPair, equalStringIteratorPair> SwarmingInvertedIndex;
        typedef RollingIndices<SwarmingInvertedIndex> SwarmingIndices;

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
        void add_candidate_counts(const numSeqs_t ampl_id, lenSeqs_t partner_len, std::vector<numSeqs_t>& cand_cnts);

        /*
         * Verify the candidates by searching for amplicons with sufficient segment matches,
         * applying the segment filter in the other direction as well,
         * and then computing the distance using the distance function on the remaining candidates.
         */
        void verify_candidates(const numSeqs_t ampl_id, const lenSeqs_t partner_len, std::vector<numSeqs_t>& cand_cnts,
                std::vector<Partner>& partners);


        const AmpliconCollection& ac_; // amplicon collection which is supported
        Distance* dist_fun_; // distance function used to find partners
        std::set<std::string> different_seqs_; // recorded different sequences

        lenSeqs_t threshold_; // threshold applied to partner search
        lenSeqs_t num_extra_segments_; // number of extra segments for segment filter
        bool break_swarms_; // flag indicating whether swarm breaking should be applied

        SwarmingIndices indices_; // inverted indices of the segment filter
        std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>> substrs_archive_; // substrings archive for more efficient segment filter

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
        typedef LabelLinkBinaryRelation<StringIteratorPair, numSeqs_t, hashStringIteratorPair, equalStringIteratorPair> SwarmingInvertedIndex;
        typedef RollingIndices<SwarmingInvertedIndex> SwarmingIndices;

    public:
        /*
         * Prepare the auxiliary data for the clustering phase (inverted indices on all amplicons).
         */
        ScoreSegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage, const numSeqs_t pool_id,
                lenSeqs_t threshold, lenSeqs_t num_extra_segments, const Configuration& config);

        /*
         * Prepare the auxiliary data for the refinement phase (inverted indices on amplicon from light swarms).
         */
        ScoreSegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                const numSeqs_t pool_id, lenSeqs_t threshold, lenSeqs_t num_extra_segments, const Configuration& config);

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
        void add_candidate_counts(const numSeqs_t ampl_id, lenSeqs_t partner_len, std::vector<numSeqs_t>& cand_cnts);

        /*
         * Verify the candidates by searching for amplicons with sufficient segment matches
         * and then computing the distance using the distance function.
         */
        void verify_candidates(const numSeqs_t ampl_id, std::vector<numSeqs_t>& cand_cnts, std::vector<Partner>& partners);


        const AmpliconCollection& ac_; // amplicon collection which is supported
        Distance* dist_fun_; // distance function used to find partners
        std::set<std::string> different_seqs_; // recorded different sequences

        lenSeqs_t threshold_; // threshold applied to partner search
        lenSeqs_t num_threshold_segments_; // number of segments based on threshold and minimum penalty for edit operation for segment filter
        lenSeqs_t num_extra_segments_; // number of extra segments for segment filter
        bool break_swarms_; // flag indicating whether swarm breaking should be applied

        SwarmingIndices indices_; // inverted indices of the segment filter
        std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>> substrs_archive_; // substrings archive for more efficient segment filter

    };

}

#endif //GEFAST_AUXILIARYDATA_HPP
