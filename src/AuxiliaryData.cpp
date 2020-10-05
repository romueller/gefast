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

#include "../include/AuxiliaryData.hpp"

namespace GeFaST {

    /* === SimpleAuxiliaryData === */

    SimpleAuxiliaryData::SimpleAuxiliaryData(const AmpliconStorage& amplicon_storage, const numSeqs_t pool_id,
            const dist_t threshold, const Configuration& config) :
            ac_(amplicon_storage.get_pool(pool_id)) {

        for (numSeqs_t i = 0; i < amplicon_storage.get_pool(pool_id).size(); i++) {
            unswarmed_.insert(i);
        }

        dist_fun_ = config.build_distance_function(amplicon_storage, threshold);
        threshold_ = threshold;
        break_swarms_ = config.break_swarms;

    }

    SimpleAuxiliaryData::SimpleAuxiliaryData(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
            const numSeqs_t pool_id, const dist_t threshold, const Configuration& config) :
            ac_(amplicon_storage.get_pool(pool_id)) {

        auto& swarms = swarm_storage.get_swarms(pool_id);

        // index all amplicons of light swarms
        for (numSeqs_t s = 0; s < swarms.size(); s++) {

            auto& swarm = swarms.get(s);

            if (swarm.mass() < config.boundary) {

                for (numSeqs_t i = 0; i < swarm.size(); i++) {
                    unswarmed_.insert(swarm.member(i));
                }

            }

        }

        dist_fun_ = config.build_distance_function(amplicon_storage, threshold);
        threshold_ = threshold;
        break_swarms_ = false;

    }

    SimpleAuxiliaryData::~SimpleAuxiliaryData() {
        delete dist_fun_;
    }

    SimpleAuxiliaryData::SimpleAuxiliaryData(const SimpleAuxiliaryData& other) : // copy constructor
            ac_(other.ac_), dist_fun_(other.dist_fun_->clone()), different_seqs_(other.different_seqs_),
            threshold_(other.threshold_), break_swarms_(other.break_swarms_), unswarmed_(other.unswarmed_) {

        // nothing else to do

    }

    SimpleAuxiliaryData::SimpleAuxiliaryData(SimpleAuxiliaryData&& other) noexcept : // move constructor
            ac_(other.ac_), dist_fun_(other.dist_fun_), different_seqs_(std::move(other.different_seqs_)),
            threshold_(other.threshold_), break_swarms_(other.break_swarms_), unswarmed_(std::move(other.unswarmed_)) {

        other.dist_fun_ = nullptr;

    }

    SimpleAuxiliaryData* SimpleAuxiliaryData::clone() const {
        return new SimpleAuxiliaryData(*this);
    }

    void SimpleAuxiliaryData::tick_off_amplicon(numSeqs_t ampl_id) {
        unswarmed_.erase(ampl_id);
    }

    bool SimpleAuxiliaryData::record_amplicon(numSeqs_t ampl_id) {
        return different_seqs_.insert(ac_.seq(ampl_id)).second;
    }

    void SimpleAuxiliaryData::clear_amplicon_records() {
        different_seqs_.clear();
    }

    std::vector<AuxiliaryData::Partner> SimpleAuxiliaryData::find_partners(numSeqs_t ampl_id) {

        std::vector<Partner> partners;

        for (auto partner_id : unswarmed_) {

            if (!break_swarms_ || (ac_.ab(ampl_id) >= ac_.ab(partner_id))) {

                auto dist = dist_fun_->distance(ac_, ampl_id, partner_id);

                if (dist <= threshold_) {
                    partners.emplace_back(partner_id, dist);
                }

            }

        }

        return partners;

    }



    /* === Support data structures and functions for segment-filter auxiliary data === */

    Substrings::Substrings(lenSeqs_t fp, lenSeqs_t lp, lenSeqs_t l) : first(fp), last(lp), len(l) {
        // nothing else to do
    }

    // Based on a combination of:
    // Li et al. (2013), A partition-based method for string similarity joins with edit-distance constraints
    // Lin et al. (2014), Large-Scale Similarity Join with Edit-Distance Constraints
    // Huang et al. (2015), A Partition-Based Bi-directional Filtering Method for String Similarity JOINs
    Substrings select_substrs(const lenSeqs_t self_len, const lenSeqs_t partner_len, const lenSeqs_t seg_index,
            const lenSeqs_t t, const lenSeqs_t k) {

        lenSeqs_t p = 0;
        lenSeqs_t l = 0;

        lenSeqs_t i = seg_index + 1; // computations from PJ-paper use segment indices from 1 to t+1

        // by even-partitioning scheme
        // * first t + k - d segments: length floor(|s| / (t+l))
        // * last d segments: length floor(|s| / (t+k)) + 1
        lenSeqs_t d = partner_len - (partner_len / (t + k)) * (t + k);
        l = (partner_len / (t + k)) + (i > t + k - d);
        p = 1 + (i - 1) * l - (i > t + k - d) * (t + k - d);


        // multimatch-aware substring selection boundaries
        /*
        // calculations from paper (assume that positions p start at 1)
        lenSeqs_t lowerL = std::max((lenSeqs_t)1, p - (i - 1));
        lenSeqs_t upperL = std::min(s.length() - l + 1, p + (i - 1));

        lenSeqs_t lowerR = std::max((lenSeqs_t)1, p + (s.length() - len) - (t + k - i));
        lenSeqs_t upperR = std::min(s.length() - l + 1, p + (s.length() - len) + (t + k - i));

        lenSeqs_t lower = std::max(lowerL, lowerR);
        lenSeqs_t upper = std::min(upperL, upperR);
        */

        /*
        // problem: underflow of non-negative type lenSeqs_t
        lenSeqs_t lower = std::max({(lenSeqs_t)1, p - (i - 1), p + (s.length() - len) - (t + k - i)});
        */
        lenSeqs_t lower = 1;
        if (p > (i - 1)) {
            lower = p - (i - 1);
        }
        if (((p + (self_len - partner_len)) > (t + k - i)) && ((p + (self_len - partner_len) - (t + k - i)) > lower)) {
            lower = p + (self_len - partner_len) - (t + k - i);
        }

        // std::min okay
        lenSeqs_t upper = std::min({
                                           self_len - l + 1,
                                           p + (i - 1),
                                           p + (self_len - partner_len) + (t + k - i)
                                   });

        return Substrings(lower - 1, upper - 1, l); // -1 for translation to index starting at 0

    }

    Substrings select_substrs_backward(const lenSeqs_t self_len, const lenSeqs_t partner_len, const lenSeqs_t seg_index,
            const lenSeqs_t t, const lenSeqs_t k) {

        lenSeqs_t p = 0;
        lenSeqs_t l = 0;

        lenSeqs_t i = seg_index + 1; // computations from PJ-paper use segment indices from 1 to t+1

        // by even-partitioning scheme
        // * first t + k - d segments: length floor(|s| / (t+l))
        // * last d segments: length floor(|s| / (t+k)) + 1
        lenSeqs_t d = partner_len - (partner_len / (t + k)) * (t + k);
        l = (partner_len / (t + k)) + (i > t + k - d);
        p = 1 + (i - 1) * l - (i > t + k - d) * (t + k - d);


        // multimatch-aware substring selection boundaries
        /*
        // calculations from paper (assume that positions p start at 1)
        lenSeqs_t lowerL = std::max((lenSeqs_t)1, p - (i - 1));
        lenSeqs_t upperL = std::min(s.length() - l + 1, p + (i - 1));

        lenSeqs_t lowerR = std::max((lenSeqs_t)1, p - (len - s.length()) - (t + k - i));
        lenSeqs_t upperR = std::min(s.length() - l + 1, p + (len - s.length()) + (t + k - i));

        lenSeqs_t lower = std::max(lowerL, lowerR);
        lenSeqs_t upper = std::min(upperL, upperR);
        */

        /*
        // problem: underflow of non-negative type lenSeqs_t
        lenSeqs_t lower = std::max({(lenSeqs_t)1, p - (i - 1), p + (s.length() - len) - (t + k - i)});
        */
        lenSeqs_t lower = 1;
        if (p > (i - 1)) {
            lower = p - (i - 1);
        }
        if ((p > ((partner_len - self_len) + (t + k - i))) && ((p - (partner_len - self_len) - (t + k - i)) > lower)) {
            lower = p - (partner_len - self_len) - (t + k - i);
        }

        // std::min okay
        lenSeqs_t upper = std::min({
                                           self_len - l + 1,
                                           p + (i - 1),
                                           p - (partner_len - self_len) + (t + k - i)
                                   });

        return Substrings(lower - 1, upper - 1, l); // -1 for translation to index starting at 0

    }

    // The first t + k - d segments have length floor(seqLen / (t + k)), while the last d segments have length ceil(seqLen / (t + k)).
    // Since d is calculated by seqLen - floor(seqLen / (t + k)) * (t + k), longer segments exist only if seqLen is not divisible by (t + k)
    // and their lengths then higher by exactly one.
    void select_segments(Segments& segments, const lenSeqs_t seq_len, const lenSeqs_t t, const lenSeqs_t k) {

        lenSeqs_t d = seq_len - (seq_len / (t + k)) * (t + k);
        lenSeqs_t p = 0;
        lenSeqs_t i = 0;

        lenSeqs_t l = (seq_len / (t + k));
        for (; i < t + k - d; i++) {

            segments[i] = std::make_pair(p, l);
            p += l;

        }

        l++;
        for (; i < t + k; i++) {

            segments[i] = std::make_pair(p, l);
            p += l;

        }

    }
    
    //TODO? custom hash function (e.g. FNV) to avoid construction of temporary string
    size_t hashStringIteratorPair::operator()(const StringIteratorPair& p) const {
        return std::hash<std::string>{}(std::string(p.first, p.second));
    }

    bool equalStringIteratorPair::operator()(const StringIteratorPair& lhs, const StringIteratorPair& rhs) const {
        return ((lhs.second - lhs.first) == (rhs.second - rhs.first)) && std::equal(lhs.first, lhs.second, rhs.first);
    }

    bool lessStringIteratorPair::operator()(const StringIteratorPair& a, const StringIteratorPair& b) const {
        return std::lexicographical_compare(a.first, a.second, b.first, b.second);
    }


    /* === SegmentFilterAuxiliaryData === */

    SegmentFilterAuxiliaryData::SegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage, const numSeqs_t pool_id,
            lenSeqs_t threshold, lenSeqs_t num_extra_segments, const Configuration& config) :
            indices_(SwarmingIndices(threshold + 1, threshold + num_extra_segments, true, false)),
            ac_(amplicon_storage.get_pool(pool_id)) {

        threshold_ = threshold;
        num_extra_segments_ = num_extra_segments;
        break_swarms_ = config.break_swarms;

        prepare();

        dist_fun_ = config.build_distance_function(amplicon_storage, threshold);

    }

    SegmentFilterAuxiliaryData::SegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, const numSeqs_t pool_id, lenSeqs_t threshold,
            lenSeqs_t num_extra_segments, const Configuration& config) :
            indices_(SwarmingIndices(threshold + 1, threshold + num_extra_segments, true, false)),
            ac_(amplicon_storage.get_pool(pool_id)) {

        threshold_ = threshold;
        num_extra_segments_ = num_extra_segments;
        break_swarms_ = false;

        prepare_light(swarm_storage.get_swarms(pool_id), config.boundary);

        dist_fun_ = config.build_distance_function(amplicon_storage, threshold);

    }

    SegmentFilterAuxiliaryData::~SegmentFilterAuxiliaryData() {
        delete dist_fun_;
    }

    SegmentFilterAuxiliaryData::SegmentFilterAuxiliaryData(const SegmentFilterAuxiliaryData& other) : // copy constructor
            ac_(other.ac_), dist_fun_(other.dist_fun_->clone()), different_seqs_(other.different_seqs_),
            threshold_(other.threshold_), num_extra_segments_(other.num_extra_segments_),
            break_swarms_(other.break_swarms_), indices_(other.indices_), substrs_archive_(other.substrs_archive_) {

        // nothing else to do

    }

    SegmentFilterAuxiliaryData::SegmentFilterAuxiliaryData(SegmentFilterAuxiliaryData&& other) noexcept : // move constructor
            ac_(other.ac_), dist_fun_(other.dist_fun_), different_seqs_(std::move(other.different_seqs_)),
            threshold_(other.threshold_), num_extra_segments_(other.num_extra_segments_),
            break_swarms_(other.break_swarms_), indices_(std::move(other.indices_)),
            substrs_archive_(std::move(other.substrs_archive_)) {

        other.dist_fun_ = nullptr;

    }

    SegmentFilterAuxiliaryData* SegmentFilterAuxiliaryData::clone() const {
        return new SegmentFilterAuxiliaryData(*this);
    }


    void SegmentFilterAuxiliaryData::tick_off_amplicon(numSeqs_t ampl_id) {

        auto& invs = indices_.get_indices_row(ac_.len(ampl_id));
        for (auto i = 0; i < threshold_ + num_extra_segments_; i++) {
            invs[i].remove_label(ampl_id);
        }

    }

    bool SegmentFilterAuxiliaryData::record_amplicon(numSeqs_t ampl_id) {
        return different_seqs_.insert(ac_.seq(ampl_id)).second;
    }

    void SegmentFilterAuxiliaryData::clear_amplicon_records() {
        different_seqs_.clear();
    }

    std::vector<AuxiliaryData::Partner> SegmentFilterAuxiliaryData::find_partners(numSeqs_t ampl_id) {

        std::vector<Partner> partners;
        std::vector<numSeqs_t> cand_cnts;

        for (lenSeqs_t partner_len = (ac_.len(ampl_id) > threshold_) * (ac_.len(ampl_id) - threshold_);
                partner_len <= ac_.len(ampl_id) + threshold_; partner_len++) {

            cand_cnts.clear();

            add_candidate_counts(ampl_id, partner_len, cand_cnts);
            verify_candidates(ampl_id, cand_cnts, partners);

        }

        return partners;

    }

    void SegmentFilterAuxiliaryData::prepare() {

        // index all amplicons
        for (numSeqs_t ampl_id = 0; ampl_id < ac_.size(); ampl_id++) {

            lenSeqs_t seq_len = ac_.len(ampl_id);

            if (!indices_.contains(seq_len)) {

                // inverted index
                indices_.roll(seq_len);

                // substrings information
                std::unordered_map<lenSeqs_t, std::vector<Substrings>>& substrs = substrs_archive_[seq_len];
                for (lenSeqs_t partner_len = (seq_len > threshold_) * (seq_len - threshold_);
                        partner_len <= seq_len + threshold_; partner_len++) {

                    std::vector<Substrings>& vec = substrs[partner_len];
                    for (lenSeqs_t segment_index = 0; segment_index < threshold_ + num_extra_segments_; segment_index++) {
                        if (partner_len <= seq_len) {
                            vec.push_back(select_substrs(seq_len, partner_len, segment_index, threshold_, num_extra_segments_));
                        } else {
                            vec.push_back(select_substrs_backward(seq_len, partner_len, segment_index, threshold_, num_extra_segments_));
                        }
                    }

                }

            }

            lenSeqs_t d = seq_len - (seq_len / (threshold_ + num_extra_segments_)) * (threshold_ + num_extra_segments_);
            lenSeqs_t p = 0;
            lenSeqs_t l = (seq_len / (threshold_ + num_extra_segments_));

            for (lenSeqs_t s = 0; s < threshold_ + num_extra_segments_; s++) {

                indices_.get_index(seq_len, s).add(StringIteratorPair(
                        ac_.seq(ampl_id) + p,
                        ac_.seq(ampl_id) + p + l + (s >= (threshold_ + num_extra_segments_ - d))),
                                                   ampl_id);

                p += l + (s >= (threshold_ + num_extra_segments_ - d));

            }

        }

    }

    void SegmentFilterAuxiliaryData::prepare_light(const Swarms& swarms, const numSeqs_t boundary) {

        std::vector<numSeqs_t> light_amplicons;

        // prepare indexing of amplicons of light swarms
        for (numSeqs_t s = 0; s < swarms.size(); s++) {

            auto& swarm = swarms.get(s);

            if (swarm.mass() < boundary && !swarm.is_attached()) {

                for (numSeqs_t i = 0; i < swarm.size(); i++) {
                    light_amplicons.push_back(swarm.member(i));
                }

            } else {

                for (numSeqs_t i = 0; i < swarm.size(); i++) {

                    lenSeqs_t seq_len = ac_.len(swarm.member(i));

                    if (substrs_archive_.find(seq_len) == substrs_archive_.end()) {

                        // substrings information of amplicons of heavy swarms
                        std::unordered_map<lenSeqs_t, std::vector<Substrings>>& substrs = substrs_archive_[seq_len];
                        for (lenSeqs_t partner_len = (seq_len > threshold_) * (seq_len - threshold_);
                                partner_len <= seq_len + threshold_; partner_len++) {

                            std::vector<Substrings>& vec = substrs[partner_len];
                            for (lenSeqs_t segment_index = 0; segment_index < threshold_ + num_extra_segments_; segment_index++) {
                                if (partner_len <= seq_len) {
                                    vec.push_back(select_substrs(seq_len, partner_len, segment_index, threshold_, num_extra_segments_));
                                } else {
                                    vec.push_back(select_substrs_backward(seq_len, partner_len, segment_index, threshold_, num_extra_segments_));
                                }
                            }

                        }

                    }

                }

            }

        }

        // sort "labels" (pool-internal integer ids) of amplicons from light swarms
        std::sort(light_amplicons.begin(), light_amplicons.end());

        // index amplicons from light swarms into relation
        for (numSeqs_t ampl_id : light_amplicons) {

            lenSeqs_t seq_len = ac_.len(ampl_id);

            if (!indices_.contains(seq_len)) { // new inverted index row
                indices_.roll(seq_len);
            }

            lenSeqs_t d = seq_len - (seq_len / (threshold_ + num_extra_segments_)) * (threshold_ + num_extra_segments_);
            lenSeqs_t p = 0;
            lenSeqs_t l = (seq_len / (threshold_ + num_extra_segments_));

            for (lenSeqs_t seg = 0; seg < threshold_ + num_extra_segments_; seg++) {

                indices_.get_index(seq_len, seg).add(StringIteratorPair(
                        ac_.seq(ampl_id) + p,
                        ac_.seq(ampl_id) + p + l + (seg >= (threshold_ + num_extra_segments_ - d))),
                                                     ampl_id);

                p += l + (seg >= (threshold_ + num_extra_segments_ - d));

            }

        }

    }

    void SegmentFilterAuxiliaryData::add_candidate_counts(const numSeqs_t ampl_id, lenSeqs_t partner_len, std::vector<numSeqs_t>& cand_cnts) {

        StringIteratorPair sip;
        for (lenSeqs_t i = 0; i < threshold_ + num_extra_segments_; i++) {

            const Substrings& subs = substrs_archive_[ac_.len(ampl_id)][partner_len][i];
            SwarmingInvertedIndex& inv = indices_.get_index(partner_len, i);
            sip.first = ac_.seq(ampl_id) + subs.first;
            sip.second = sip.first + subs.len;

            for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++, sip.first++, sip.second++) {
                inv.add_label_counts_of(sip, cand_cnts);
            }

        }

    }

    void SegmentFilterAuxiliaryData::verify_candidates(const numSeqs_t ampl_id, std::vector<numSeqs_t>& cand_cnts,
            std::vector<Partner>& partners) {

        std::sort(cand_cnts.begin(), cand_cnts.end());

        // general pigeonhole principle: for being a candidate, at least k segments have to be matched
        lenSeqs_t cnt = 0;
        numSeqs_t prev_cand = (cand_cnts.empty()) ? 0 : cand_cnts.front();
        for (auto cand_id : cand_cnts) {

            if (prev_cand != cand_id) {

                if ((cnt >= num_extra_segments_) && (!break_swarms_ || ac_.ab(ampl_id) >= ac_.ab(prev_cand))) {

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

        if ((cnt >= num_extra_segments_) && (!break_swarms_ || ac_.ab(ampl_id) >= ac_.ab(prev_cand))) {

            auto dist = dist_fun_->distance(ac_, ampl_id, prev_cand);

            if (dist <= threshold_) {
                partners.emplace_back(prev_cand, dist);
            }

        }

    }


    /* === TwoWaySegmentFilterAuxiliaryData === */

    TwoWaySegmentFilterAuxiliaryData::TwoWaySegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage,
            const numSeqs_t pool_id, lenSeqs_t threshold, lenSeqs_t num_extra_segments, const Configuration& config) :
            indices_(SwarmingIndices(threshold + 1, threshold + num_extra_segments, true, false)),
            ac_(amplicon_storage.get_pool(pool_id)) {

        threshold_ = threshold;
        num_extra_segments_ = num_extra_segments;
        break_swarms_ = config.break_swarms;

        prepare();

        dist_fun_ = config.build_distance_function(amplicon_storage, threshold);

        segments_ = Segments(threshold_ + num_extra_segments_);
        segment_strs_ = std::vector<std::string>(threshold_ + num_extra_segments_);

    }

    TwoWaySegmentFilterAuxiliaryData::TwoWaySegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, const numSeqs_t pool_id, lenSeqs_t threshold,
            lenSeqs_t num_extra_segments, const Configuration& config) :
            indices_(SwarmingIndices(threshold + 1, threshold + num_extra_segments, true, false)),
            ac_(amplicon_storage.get_pool(pool_id)) {

        threshold_ = threshold;
        num_extra_segments_ = num_extra_segments;
        break_swarms_ = false;

        prepare_light(swarm_storage.get_swarms(pool_id), config.boundary);

        dist_fun_ = config.build_distance_function(amplicon_storage, threshold);

        segments_ = Segments(threshold_ + num_extra_segments_);
        segment_strs_ = std::vector<std::string>(threshold_ + num_extra_segments_);

    }

    TwoWaySegmentFilterAuxiliaryData::~TwoWaySegmentFilterAuxiliaryData() {
        delete dist_fun_;
    }

    TwoWaySegmentFilterAuxiliaryData::TwoWaySegmentFilterAuxiliaryData(const TwoWaySegmentFilterAuxiliaryData& other) : // copy constructor
            ac_(other.ac_), dist_fun_(other.dist_fun_->clone()), different_seqs_(other.different_seqs_),
            threshold_(other.threshold_), num_extra_segments_(other.num_extra_segments_),
            break_swarms_(other.break_swarms_), indices_(other.indices_), substrs_archive_(other.substrs_archive_),
            segments_(other.segments_), segment_strs_(other.segment_strs_) {

        // nothing else to do

    }

    TwoWaySegmentFilterAuxiliaryData::TwoWaySegmentFilterAuxiliaryData(TwoWaySegmentFilterAuxiliaryData&& other) noexcept : // move constructor
            ac_(other.ac_), dist_fun_(other.dist_fun_), different_seqs_(std::move(other.different_seqs_)),
            threshold_(other.threshold_), num_extra_segments_(other.num_extra_segments_),
            break_swarms_(other.break_swarms_), indices_(std::move(other.indices_)),
            substrs_archive_(std::move(other.substrs_archive_)), segments_(std::move(other.segments_)),
            segment_strs_(std::move(other.segment_strs_)) {

        other.dist_fun_ = nullptr;

    }

    TwoWaySegmentFilterAuxiliaryData* TwoWaySegmentFilterAuxiliaryData::clone() const {
        return new TwoWaySegmentFilterAuxiliaryData(*this);
    }

    void TwoWaySegmentFilterAuxiliaryData::tick_off_amplicon(numSeqs_t ampl_id) {

        auto& invs = indices_.get_indices_row(ac_.len(ampl_id));
        for (auto i = 0; i < threshold_ + num_extra_segments_; i++) {
            invs[i].remove_label(ampl_id);
        }

    }

    bool TwoWaySegmentFilterAuxiliaryData::record_amplicon(numSeqs_t ampl_id) {
        return different_seqs_.insert(ac_.seq(ampl_id)).second;
    }

    void TwoWaySegmentFilterAuxiliaryData::clear_amplicon_records() {
        different_seqs_.clear();
    }

    std::vector<AuxiliaryData::Partner> TwoWaySegmentFilterAuxiliaryData::find_partners(numSeqs_t ampl_id) {

        std::vector<Partner> partners;
        std::vector<numSeqs_t> cand_cnts;

        select_segments(segments_, ac_.len(ampl_id), threshold_, num_extra_segments_);
        for (auto i = 0; i < threshold_ + num_extra_segments_; i++) {
            segment_strs_[i] = std::string(ac_.seq(ampl_id) + segments_[i].first,
                    ac_.seq(ampl_id) + segments_[i].first + segments_[i].second);
        }

        for (lenSeqs_t partner_len = (ac_.len(ampl_id) > threshold_) * (ac_.len(ampl_id) - threshold_);
                partner_len <= ac_.len(ampl_id) + threshold_; partner_len++) {

            cand_cnts.clear();

            add_candidate_counts(ampl_id, partner_len, cand_cnts);
            verify_candidates(ampl_id, partner_len, cand_cnts, partners);

        }

        return partners;

    }

    void TwoWaySegmentFilterAuxiliaryData::prepare() {

        // index all amplicons
        for (numSeqs_t ampl_id = 0; ampl_id < ac_.size(); ampl_id++) {

            lenSeqs_t seq_len = ac_.len(ampl_id);

            if (!indices_.contains(seq_len)) {

                // inverted index
                indices_.roll(seq_len);

                // substrings information
                std::unordered_map<lenSeqs_t, std::vector<Substrings>>& substrs = substrs_archive_[seq_len];
                for (lenSeqs_t partner_len = (seq_len > threshold_) * (seq_len - threshold_);
                        partner_len <= seq_len + threshold_; partner_len++) {

                    std::vector<Substrings>& vec = substrs[partner_len];
                    for (lenSeqs_t segment_index = 0; segment_index < threshold_ + num_extra_segments_; segment_index++) {
                        if (partner_len <= seq_len) {
                            vec.push_back(select_substrs(seq_len, partner_len, segment_index, threshold_, num_extra_segments_));
                        } else {
                            vec.push_back(select_substrs_backward(seq_len, partner_len, segment_index, threshold_, num_extra_segments_));
                        }
                    }

                }

            }

            lenSeqs_t d = seq_len - (seq_len / (threshold_ + num_extra_segments_)) * (threshold_ + num_extra_segments_);
            lenSeqs_t p = 0;
            lenSeqs_t l = (seq_len / (threshold_ + num_extra_segments_));

            for (lenSeqs_t s = 0; s < threshold_ + num_extra_segments_; s++) {

                indices_.get_index(seq_len, s).add(StringIteratorPair(
                        ac_.seq(ampl_id) + p,
                        ac_.seq(ampl_id) + p + l + (s >= (threshold_ + num_extra_segments_ - d))),
                                                   ampl_id);

                p += l + (s >= (threshold_ + num_extra_segments_ - d));

            }

        }

    }

    void TwoWaySegmentFilterAuxiliaryData::prepare_light(const Swarms& swarms, const numSeqs_t boundary) {

        std::vector<numSeqs_t> light_amplicons;

        // prepare indexing of amplicons of light swarms
        for (numSeqs_t s = 0; s < swarms.size(); s++) {

            auto& swarm = swarms.get(s);

            if (swarm.mass() < boundary && !swarm.is_attached()) {

                for (numSeqs_t i = 0; i < swarm.size(); i++) {
                    light_amplicons.push_back(swarm.member(i));
                }

            }

            for (numSeqs_t i = 0; i < swarm.size(); i++) {

                lenSeqs_t seq_len = ac_.len(swarm.member(i));

                if (substrs_archive_.find(seq_len) == substrs_archive_.end()) {

                    // substrings information of amplicons of heavy swarms
                    std::unordered_map<lenSeqs_t, std::vector<Substrings>>& substrs = substrs_archive_[seq_len];
                    for (lenSeqs_t partner_len = (seq_len > threshold_) * (seq_len - threshold_);
                            partner_len <= seq_len + threshold_; partner_len++) {

                        std::vector<Substrings>& vec = substrs[partner_len];
                        for (lenSeqs_t segment_index = 0; segment_index < threshold_ + num_extra_segments_; segment_index++) {
                            if (partner_len <= seq_len) {
                                vec.push_back(select_substrs(seq_len, partner_len, segment_index, threshold_, num_extra_segments_));
                            } else {
                                vec.push_back(select_substrs_backward(seq_len, partner_len, segment_index, threshold_, num_extra_segments_));
                            }
                        }

                    }

                }

            }

        }

        // sort "labels" (pool-internal integer ids) of amplicons from light swarms
        std::sort(light_amplicons.begin(), light_amplicons.end());

        // index amplicons from light swarms into relation
        for (numSeqs_t ampl_id : light_amplicons) {

            lenSeqs_t seq_len = ac_.len(ampl_id);

            if (!indices_.contains(seq_len)) { // new inverted index row
                indices_.roll(seq_len);
            }

            lenSeqs_t d = seq_len - (seq_len / (threshold_ + num_extra_segments_)) * (threshold_ + num_extra_segments_);
            lenSeqs_t p = 0;
            lenSeqs_t l = (seq_len / (threshold_ + num_extra_segments_));

            for (lenSeqs_t seg = 0; seg < threshold_ + num_extra_segments_; seg++) {

                indices_.get_index(seq_len, seg).add(StringIteratorPair(
                        ac_.seq(ampl_id) + p,
                        ac_.seq(ampl_id) + p + l + (seg >= (threshold_ + num_extra_segments_ - d))),
                                                     ampl_id);

                p += l + (seg >= (threshold_ + num_extra_segments_ - d));

            }

        }

    }

    void TwoWaySegmentFilterAuxiliaryData::add_candidate_counts(const numSeqs_t ampl_id, lenSeqs_t partner_len,
            std::vector<numSeqs_t>& cand_cnts) {

        StringIteratorPair sip;
        for (lenSeqs_t i = 0; i < threshold_ + num_extra_segments_; i++) {

            const Substrings& subs = substrs_archive_[ac_.len(ampl_id)][partner_len][i];
            SwarmingInvertedIndex& inv = indices_.get_index(partner_len, i);
            sip.first = ac_.seq(ampl_id) + subs.first;
            sip.second = sip.first + subs.len;

            for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++, sip.first++, sip.second++) {
                inv.add_label_counts_of(sip, cand_cnts);
            }

        }

    }

    void TwoWaySegmentFilterAuxiliaryData::verify_candidates(const numSeqs_t ampl_id, const lenSeqs_t partner_len,
            std::vector<numSeqs_t>& cand_cnts, std::vector<Partner>& partners) {

        std::sort(cand_cnts.begin(), cand_cnts.end());

        // general pigeonhole principle: for being a candidate, at least k segments have to be matched
        std::vector<Substrings>& cand_substrs = substrs_archive_[partner_len][ac_.len(ampl_id)];

        lenSeqs_t cnt = 0;
        numSeqs_t prev_cand = (cand_cnts.empty()) ? 0 : cand_cnts.front();
        for (auto cand_id : cand_cnts) {

            if (prev_cand != cand_id) {

                if ((cnt >= num_extra_segments_) && (!break_swarms_ || ac_.ab(ampl_id) >= ac_.ab(prev_cand))) {

                    cnt = 0; // reset to 0 and count matching segments in the other direction
                    std::string cand_str = std::string(ac_.seq(prev_cand));
                    for (lenSeqs_t i = 0; i < threshold_ + num_extra_segments_ && cnt < num_extra_segments_; i++) {
                        cnt += (cand_str.substr(
                                cand_substrs[i].first,
                                (cand_substrs[i].last - cand_substrs[i].first) + cand_substrs[i].len
                        ).find(segment_strs_[i]) < std::string::npos);
                    }

                    if (cnt == num_extra_segments_) {

                        auto dist = dist_fun_->distance(ac_, ampl_id, prev_cand);

                        if (dist <= threshold_) {
                            partners.emplace_back(prev_cand, dist);
                        }

                    }

                }

                cnt = 1;
                prev_cand = cand_id;


            } else {
                cnt++;
            }

        }

        if ((cnt >= num_extra_segments_) && (!break_swarms_ || ac_.ab(ampl_id) >= ac_.ab(prev_cand))) {

            cnt = 0; // reset to 0 and count matching segments in the other direction
            std::string cand_str = std::string(ac_.seq(prev_cand));
            for (lenSeqs_t i = 0; i < threshold_ + num_extra_segments_ && cnt < num_extra_segments_; i++) {
                cnt += (cand_str.substr(
                        cand_substrs[i].first,
                        (cand_substrs[i].last - cand_substrs[i].first) + cand_substrs[i].len
                ).find(segment_strs_[i]) < std::string::npos);
            }

            if (cnt == num_extra_segments_) {

                auto dist = dist_fun_->distance(ac_, ampl_id, prev_cand);

                if (dist <= threshold_) {
                    partners.emplace_back(prev_cand, dist);
                }

            }

        }

    }


    /* === ScoreSegmentFilterAuxiliaryData === */

    ScoreSegmentFilterAuxiliaryData::ScoreSegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage,
            const numSeqs_t pool_id, lenSeqs_t threshold, lenSeqs_t num_extra_segments, const Configuration& config) :
            indices_(SwarmingIndices(threshold + 1, threshold + num_extra_segments, true, false)),
            ac_(amplicon_storage.get_pool(pool_id)) {

        ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);

        threshold_ = threshold;
        num_threshold_segments_ = threshold / std::min({sf.penalty_mismatch, sf.penalty_open, sf.penalty_extend});
        num_extra_segments_ = num_extra_segments;
        break_swarms_ = config.break_swarms;

        prepare();

        dist_fun_ = config.build_distance_function(amplicon_storage, threshold);

    }

    ScoreSegmentFilterAuxiliaryData::ScoreSegmentFilterAuxiliaryData(const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, const numSeqs_t pool_id, lenSeqs_t threshold,
            lenSeqs_t num_extra_segments, const Configuration& config) :
            indices_(SwarmingIndices(threshold + 1, threshold + num_extra_segments, true, false)),
            ac_(amplicon_storage.get_pool(pool_id)) {

        ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);

        threshold_ = threshold;
        num_threshold_segments_ = threshold / std::min({sf.penalty_mismatch, sf.penalty_open, sf.penalty_extend});
        num_extra_segments_ = num_extra_segments;
        break_swarms_ = false;

        prepare_light(swarm_storage.get_swarms(pool_id), config.boundary);

        dist_fun_ = config.build_distance_function(amplicon_storage, threshold);

    }

    ScoreSegmentFilterAuxiliaryData::~ScoreSegmentFilterAuxiliaryData() {
        delete dist_fun_;
    }

    ScoreSegmentFilterAuxiliaryData::ScoreSegmentFilterAuxiliaryData(const ScoreSegmentFilterAuxiliaryData& other) : // copy constructor
        ac_(other.ac_), dist_fun_(other.dist_fun_->clone()), different_seqs_(other.different_seqs_),
        threshold_(other.threshold_), num_threshold_segments_(other.num_threshold_segments_),
        num_extra_segments_(other.num_extra_segments_), break_swarms_(other.break_swarms_),
        indices_(other.indices_), substrs_archive_(other.substrs_archive_) {

        // nothing else to do

    }

    ScoreSegmentFilterAuxiliaryData::ScoreSegmentFilterAuxiliaryData(ScoreSegmentFilterAuxiliaryData&& other) noexcept : // move constructor
        ac_(other.ac_), dist_fun_(other.dist_fun_), different_seqs_(std::move(other.different_seqs_)),
        threshold_(other.threshold_), num_threshold_segments_(other.num_threshold_segments_),
        num_extra_segments_(other.num_extra_segments_), break_swarms_(other.break_swarms_),
        indices_(std::move(other.indices_)), substrs_archive_(std::move(other.substrs_archive_)) {

        other.dist_fun_ = nullptr;

    }

    ScoreSegmentFilterAuxiliaryData* ScoreSegmentFilterAuxiliaryData::clone() const {
        return new ScoreSegmentFilterAuxiliaryData(*this);
    }

    void ScoreSegmentFilterAuxiliaryData::tick_off_amplicon(numSeqs_t ampl_id) {

        auto& invs = indices_.get_indices_row(ac_.len(ampl_id));
        for (auto i = 0; i < num_threshold_segments_ + num_extra_segments_; i++) {
            invs[i].remove_label(ampl_id);
        }

    }

    bool ScoreSegmentFilterAuxiliaryData::record_amplicon(numSeqs_t ampl_id) {
        return different_seqs_.insert(ac_.seq(ampl_id)).second;
    }

    void ScoreSegmentFilterAuxiliaryData::clear_amplicon_records() {
        different_seqs_.clear();
    }

    std::vector<AuxiliaryData::Partner> ScoreSegmentFilterAuxiliaryData::find_partners(numSeqs_t ampl_id) {

        std::vector<Partner> partners;
        std::vector<numSeqs_t> cand_cnts;

        for (lenSeqs_t partner_len = (ac_.len(ampl_id) > num_threshold_segments_) * (ac_.len(ampl_id) - num_threshold_segments_);
                partner_len <= ac_.len(ampl_id) + num_threshold_segments_; partner_len++) {

            cand_cnts.clear();

            add_candidate_counts(ampl_id, partner_len, cand_cnts);
            verify_candidates(ampl_id, cand_cnts, partners);

        }

        return partners;

    }

    void ScoreSegmentFilterAuxiliaryData::prepare() {

        // index all amplicons
        for (numSeqs_t ampl_id = 0; ampl_id < ac_.size(); ampl_id++) {

            lenSeqs_t seq_len = ac_.len(ampl_id);

            if (!indices_.contains(seq_len)) {

                // inverted index
                indices_.roll(seq_len);

                // substrings information
                std::unordered_map<lenSeqs_t, std::vector<Substrings>>& substrs = substrs_archive_[seq_len];
                for (lenSeqs_t partner_len = (seq_len > num_threshold_segments_) * (seq_len - num_threshold_segments_);
                        partner_len <= seq_len + num_threshold_segments_; partner_len++) {

                    std::vector<Substrings>& vec = substrs[partner_len];
                    for (lenSeqs_t segment_index = 0; segment_index < num_threshold_segments_ + num_extra_segments_; segment_index++) {
                        if (partner_len <= seq_len) {
                            vec.push_back(select_substrs(seq_len, partner_len, segment_index, num_threshold_segments_, num_extra_segments_));
                        } else {
                            vec.push_back(select_substrs_backward(seq_len, partner_len, segment_index, num_threshold_segments_, num_extra_segments_));
                        }
                    }

                }

            }

            lenSeqs_t d = seq_len - (seq_len / (num_threshold_segments_ + num_extra_segments_)) * (num_threshold_segments_ + num_extra_segments_);
            lenSeqs_t p = 0;
            lenSeqs_t l = (seq_len / (num_threshold_segments_ + num_extra_segments_));

            for (lenSeqs_t s = 0; s < num_threshold_segments_ + num_extra_segments_; s++) {

                indices_.get_index(seq_len, s).add(StringIteratorPair(
                        ac_.seq(ampl_id) + p,
                        ac_.seq(ampl_id) + p + l + (s >= (num_threshold_segments_ + num_extra_segments_ - d))),
                                                   ampl_id);

                p += l + (s >= (num_threshold_segments_ + num_extra_segments_ - d));

            }

        }

    }

    void ScoreSegmentFilterAuxiliaryData::prepare_light(const Swarms& swarms, const numSeqs_t boundary) {

        std::vector<numSeqs_t> light_amplicons;

        // prepare indexing of amplicons of light swarms
        for (numSeqs_t s = 0; s < swarms.size(); s++) {

            auto& swarm = swarms.get(s);

            if (swarm.mass() < boundary && !swarm.is_attached()) {

                for (numSeqs_t i = 0; i < swarm.size(); i++) {
                    light_amplicons.push_back(swarm.member(i));
                }

            } else {

                for (numSeqs_t i = 0; i < swarm.size(); i++) {

                    lenSeqs_t seq_len = ac_.len(swarm.member(i));

                    if (substrs_archive_.find(seq_len) == substrs_archive_.end()) {

                        // substrings information of amplicons of heavy swarms
                        std::unordered_map<lenSeqs_t, std::vector<Substrings>>& substrs = substrs_archive_[seq_len];
                        for (lenSeqs_t partner_len = (seq_len > num_threshold_segments_) * (seq_len - num_threshold_segments_);
                                partner_len <= seq_len + num_threshold_segments_; partner_len++) {

                            std::vector<Substrings>& vec = substrs[partner_len];
                            for (lenSeqs_t segment_index = 0; segment_index < num_threshold_segments_ + num_extra_segments_; segment_index++) {
                                if (partner_len <= seq_len) {
                                    vec.push_back(select_substrs(seq_len, partner_len, segment_index, num_threshold_segments_, num_extra_segments_));
                                } else {
                                    vec.push_back(select_substrs_backward(seq_len, partner_len, segment_index, num_threshold_segments_, num_extra_segments_));
                                }
                            }

                        }

                    }

                }

            }

        }

        // sort "labels" (pool-internal integer ids) of amplicons from light swarms
        std::sort(light_amplicons.begin(), light_amplicons.end());

        // index amplicons from light swarms into relation
        for (numSeqs_t ampl_id : light_amplicons) {

            lenSeqs_t seq_len = ac_.len(ampl_id);

            if (!indices_.contains(seq_len)) { // new inverted index row
                indices_.roll(seq_len);
            }

            lenSeqs_t d = seq_len - (seq_len / (num_threshold_segments_ + num_extra_segments_)) * (num_threshold_segments_ + num_extra_segments_);
            lenSeqs_t p = 0;
            lenSeqs_t l = (seq_len / (num_threshold_segments_ + num_extra_segments_));

            for (lenSeqs_t seg = 0; seg < num_threshold_segments_ + num_extra_segments_; seg++) {

                indices_.get_index(seq_len, seg).add(StringIteratorPair(
                        ac_.seq(ampl_id) + p,
                        ac_.seq(ampl_id) + p + l + (seg >= (num_threshold_segments_ + num_extra_segments_ - d))),
                                                     ampl_id);

                p += l + (seg >= (num_threshold_segments_ + num_extra_segments_ - d));

            }

        }

    }

    void ScoreSegmentFilterAuxiliaryData::add_candidate_counts(const numSeqs_t ampl_id, lenSeqs_t partner_len,
            std::vector<numSeqs_t>& cand_cnts) {

        StringIteratorPair sip;
        for (lenSeqs_t i = 0; i < num_threshold_segments_ + num_extra_segments_; i++) {

            const Substrings& subs = substrs_archive_[ac_.len(ampl_id)][partner_len][i];
            SwarmingInvertedIndex& inv = indices_.get_index(partner_len, i);
            sip.first = ac_.seq(ampl_id) + subs.first;
            sip.second = sip.first + subs.len;

            for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++, sip.first++, sip.second++) {
                inv.add_label_counts_of(sip, cand_cnts);
            }

        }

    }

    void ScoreSegmentFilterAuxiliaryData::verify_candidates(const numSeqs_t ampl_id, std::vector<numSeqs_t>& cand_cnts,
            std::vector<Partner>& partners) {

        std::sort(cand_cnts.begin(), cand_cnts.end());

        // general pigeonhole principle: for being a candidate, at least k segments have to be matched
        lenSeqs_t cnt = 0;
        numSeqs_t prev_cand = (cand_cnts.empty()) ? 0 : cand_cnts.front();
        for (auto cand_id : cand_cnts) {

            if (prev_cand != cand_id) {

                if ((cnt >= num_extra_segments_) && (!break_swarms_ || ac_.ab(ampl_id) >= ac_.ab(prev_cand))) {

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

        if ((cnt >= num_extra_segments_) && (!break_swarms_ || ac_.ab(ampl_id) >= ac_.ab(prev_cand))) {

            auto dist = dist_fun_->distance(ac_, ampl_id, prev_cand);

            if (dist <= threshold_) {
                partners.emplace_back(prev_cand, dist);
            }

        }

    }

}