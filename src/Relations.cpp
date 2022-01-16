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

#include "../include/Relations.hpp"

namespace GeFaST {

    SegmentRelations::SegmentRelations(const AmpliconCollection& ac, numSeqs_t num_segs) {

        pool_size_ = ac.size();
        num_segs_ = num_segs;

        positions_ = new position_type[pool_size_ * num_segs_];

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

    }

    SegmentRelations::~SegmentRelations() {

        delete[] access_sip_;
        delete[] positions_;

    }

    SegmentRelations::SegmentRelations(const SegmentRelations& other) : // copy constructor
            pool_size_(other.pool_size_), num_segs_(other.num_segs_), access_len_(other.access_len_) {

        positions_ = new position_type[pool_size_ * num_segs_];
        std::copy(other.positions_, other.positions_ + pool_size_ * num_segs_, positions_);

        access_sip_ = new sip_map_type[access_len_.size() * num_segs_];
        for (auto i = 0; i < access_len_.size() * num_segs_; i++) {
            for (auto& kv : other.access_sip_[i]) {
                access_sip_[i][std::unique_ptr<SubstringWrapper>(kv.first->clone())] = kv.second;
            }
        }

    }

    SegmentRelations::SegmentRelations(SegmentRelations&& other) noexcept : // move constructor
            pool_size_(other.pool_size_), num_segs_(other.num_segs_), access_len_(std::move(other.access_len_)) {

        positions_ = other.positions_; other.positions_ = nullptr;
        access_sip_ = other.access_sip_; other.access_sip_ = nullptr;

    }

    SegmentRelations& SegmentRelations::operator=(const SegmentRelations& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] access_sip_;
        delete[] positions_;

        // copy new resources
        pool_size_ = other.pool_size_;
        num_segs_ = other.num_segs_;

        positions_ = new position_type[pool_size_ * num_segs_];
        std::copy(other.positions_, other.positions_ + pool_size_ * num_segs_, positions_);

        access_len_ = other.access_len_;

        access_sip_ = new sip_map_type[access_len_.size() * num_segs_];
        for (auto i = 0; i < access_len_.size() * num_segs_; i++) {
            for (auto& kv : other.access_sip_[i]) {
                access_sip_[i][std::unique_ptr<SubstringWrapper>(kv.first->clone())] = kv.second;
            }
        }

        return *this;

    }

    SegmentRelations& SegmentRelations::operator=(SegmentRelations&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] access_sip_;
        delete[] positions_;

        // copy new resources
        pool_size_ = other.pool_size_;
        num_segs_ = other.num_segs_;

        positions_ = other.positions_; other.positions_ = nullptr;

        access_len_ = std::move(other.access_len_);
        access_sip_ = other.access_sip_; other.access_sip_ = nullptr;

        return *this;

    }

    void SegmentRelations::add_counts(const numSeqs_t len, const numSeqs_t seg, const std::unique_ptr<SubstringWrapper>& sub_sip, std::vector<numSeqs_t>& cand_cnts) {

        auto& acc = access_sip_[access_len_[len] + seg];
        auto iter = acc.find(sub_sip);

        if (iter != acc.end()) {
            cand_cnts.insert(cand_cnts.end(), iter->second.begin(), iter->second.end());
        }

    }
    void SegmentRelations::add_counts_checked(const numSeqs_t len, const numSeqs_t seg, const std::unique_ptr<SubstringWrapper>& sub_sip, std::vector<numSeqs_t>& cand_cnts) {

        if (access_len_.find(len) != access_len_.end()) {

            auto& acc = access_sip_[access_len_[len] + seg];
            auto iter = acc.find(sub_sip);

            if (iter != acc.end()) {
                cand_cnts.insert(cand_cnts.end(), iter->second.begin(), iter->second.end());
            }

        }

    }

    void SegmentRelations::add_neighboured_counts(const numSeqs_t len, const numSeqs_t seg, std::unique_ptr<SubstringWrapper>& sub_sip,
                                numSeqs_t num_neigh, std::vector<numSeqs_t>& cand_cnts) {

        auto& acc = access_sip_[access_len_[len] + seg];

        for (auto n = 0; n <= num_neigh; n++) {

            auto iter = acc.find(sub_sip);

            if (iter != acc.end()) {
                cand_cnts.insert(cand_cnts.end(), iter->second.begin(), iter->second.end());
            }

            sub_sip->shift();

        }

    }
    void SegmentRelations::add_neighboured_counts_checked(const numSeqs_t len, const numSeqs_t seg, std::unique_ptr<SubstringWrapper>& sip,
                                        numSeqs_t num_neigh, std::vector<numSeqs_t>& cand_cnts) {

        if (access_len_.find(len) != access_len_.end()) {

            auto& acc = access_sip_[access_len_[len] + seg];

            for (auto n = 0; n <= num_neigh; n++) {

                auto iter = acc.find(sip);

                if (iter != acc.end()) {
                    cand_cnts.insert(cand_cnts.end(), iter->second.begin(), iter->second.end());
                }

                sip->shift();

            }

        }

    }

    void SegmentRelations::record_segment(const numSeqs_t len, const numSeqs_t seg, const std::unique_ptr<SubstringWrapper>& sub_sip, const numSeqs_t ampl_id) {

        auto& acc = access_sip_[access_len_[len] + seg];
        acc[std::unique_ptr<SubstringWrapper>(sub_sip->clone())].emplace_back(ampl_id);
        positions_[ampl_id * num_segs_ + seg] = &(*acc.find(sub_sip));

    }

    void SegmentRelations::remove_amplicon(const numSeqs_t ampl_id) {

        for (numSeqs_t pos = ampl_id * num_segs_, seg = 0; seg < num_segs_; seg++, pos++) {

            auto& vec = positions_[pos]->second;
            vec.erase(std::find(vec.begin(), vec.end(), ampl_id)); // (linear, part) lookup + remove (halfway through)
//                vec.erase(std::remove_if(vec.begin(), vec.end(), [ampl_id](numSeqs_t a){return a == ampl_id;}), vec.end()); // (linear, full) lookup + remove (at the end)
//                vec.erase(std::lower_bound(vec.begin(), vec.end(), ampl_id)); // (binary) lookup + remove (halfway through)

        }

    }

    bool SegmentRelations::contains_length(const numSeqs_t len) const {
        return access_len_.find(len) != access_len_.end();
    }

    size_t SegmentRelations::size_in_bytes() const {

        size_t allocated_access_sip = 0;
        for (auto i = 0; i < access_len_.size() * num_segs_; i++) {
            allocated_access_sip += sip_map_multiple_size_in_bytes(access_sip_[i]);
        }

        return sizeof(SegmentRelations) // SegmentRelations itself
            + allocated_access_sip // access_sip_
            + umap_pod_size_in_bytes(access_len_) - sizeof(length_map_type) // access_len_ (do not count size of length_map_type object twice)
            + sizeof(position_type*) * pool_size_ * num_segs_; // positions_

    }

    void SegmentRelations::show_memory() const {

        size_t allocated_access_sip = 0;
        for (auto i = 0; i < access_len_.size() * num_segs_; i++) {
            allocated_access_sip += sip_map_multiple_size_in_bytes(access_sip_[i]);
        }

        std::cout << "#  * SegmentRelations " << std::endl;
        std::cout << "#  * " << std::endl;
        std::cout << "#  * sizeof(SegmentRelations): " << sizeof(SegmentRelations) << " bytes" << std::endl;
        std::cout << "#  * sizeof(sip_map_type): " << sizeof(sip_map_type) << " bytes" << std::endl;
        std::cout << "#  * sizeof(sip_map_type*): " << sizeof(sip_map_type*) << " bytes" << std::endl;
        std::cout << "#  * sizeof(length_map_type): " << sizeof(length_map_type) << " bytes" << std::endl;
        std::cout << "#  * sizeof(position_type): " << sizeof(position_type) << " bytes" << std::endl;
        std::cout << "#  * sizeof(position_type*): " << sizeof(position_type*) << " bytes" << std::endl;
        std::cout << "#  * sizeof(numSeqs_t): " << sizeof(numSeqs_t) << " bytes" << std::endl;
        std::cout << "#  * " << std::endl;
        std::cout << "#  * Total size: " << size_in_bytes() << " bytes" << std::endl;
        std::cout << "#  * access_sip_: " << allocated_access_sip << " bytes" << std::endl;
        std::cout << "#  * access_len_: " << umap_pod_size_in_bytes(access_len_) << " bytes" << std::endl;
        std::cout << "#  * positions_: " << (sizeof(position_type*) * pool_size_ * num_segs_) << " bytes" << std::endl;

    }

}