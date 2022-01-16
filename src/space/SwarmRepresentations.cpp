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

#include <numeric>

#include "../../include/space/SwarmRepresentations.hpp"

namespace GeFaST {

    /* === SharingSplitDacsSwarms === */

    SharingSplitDacsSwarms::SharingSplitDacsSwarms(const numSeqs_t pool_size, const ChunkLengths& config)
            : all_members_(config.ids, pool_size, config.level_factor, config.extend_factor, config.opt),
              all_parents_(config.ids, pool_size, config.level_factor, config.extend_factor, config.opt),
              all_parent_dists_(config.dists, pool_size, config.level_factor, config.extend_factor, config.opt),
              all_gens_(config.gens, pool_size, config.level_factor, config.extend_factor, config.opt),
              all_rads_(config.rads, pool_size, config.level_factor, config.extend_factor, config.opt) {

        cur_gen_ = 0;
        total_num_members_ = 0;
        finalised_ = false;

    }

    Swarm& SharingSplitDacsSwarms::initialise_cluster(numSeqs_t s, numSeqs_t ab) {

        transfer_generations();

        cur_members_.push_back(s);
        cur_parents_.push_back(s);
        cur_dists_.push_back(0);
        cur_gen_ = 0;
        cur_rads_.push_back(0);

        swarms_.push_back(new ForwardingDacsSwarm(*this, ab));
        return *(swarms_.back());

    }

    SharingSplitDacsSwarms::~SharingSplitDacsSwarms() {

        for (auto i = 0; i < swarms_.size(); i++) {
            delete swarms_[i]; swarms_[i] = nullptr;
        }

    }

    SharingSplitDacsSwarms::SharingSplitDacsSwarms(const SharingSplitDacsSwarms& other) :
            all_members_(other.all_members_), all_parents_(other.all_parents_), all_parent_dists_(other.all_parent_dists_),
            all_gens_(other.all_gens_), all_rads_(other.all_rads_), next_members_(other.next_members_),
            next_parents_(other.next_parents_), next_dists_(other.next_dists_),
            next_rads_(other.next_rads_), cur_members_(other.cur_members_),
            cur_parents_(other.cur_parents_), cur_dists_(other.cur_dists_), cur_rads_(other.cur_rads_) { // copy constructor

        cur_gen_ = other.cur_gen_;
        total_num_members_ = other.total_num_members_;
        finalised_ = other.finalised_;

        swarms_.reserve(other.size());
        for (auto s : other.swarms_) {
            swarms_.push_back(s->clone());
        }

    }

    SharingSplitDacsSwarms::SharingSplitDacsSwarms(SharingSplitDacsSwarms&& other) noexcept :
            swarms_(std::move(other.swarms_)), all_members_(std::move(other.all_members_)),
            all_parents_(std::move(other.all_parents_)), all_parent_dists_(std::move(other.all_parent_dists_)),
            all_gens_(std::move(other.all_gens_)), all_rads_(std::move(other.all_rads_)),
            next_members_(std::move(other.next_members_)), next_parents_(std::move(other.next_parents_)),
            next_dists_(std::move(other.next_dists_)),
            next_rads_(std::move(other.next_rads_)), cur_parents_(std::move(other.cur_parents_)),
            cur_dists_(std::move(other.cur_dists_)), cur_rads_(std::move(other.cur_rads_)) { // move constructor

        cur_gen_ = other.cur_gen_;
        total_num_members_ = other.total_num_members_;
        finalised_ = other.finalised_;

    }

    SharingSplitDacsSwarms& SharingSplitDacsSwarms::operator=(const SharingSplitDacsSwarms& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        for (auto i = 0; i < swarms_.size(); i++) {
            delete swarms_[i]; swarms_[i] = nullptr;
        }
        swarms_.clear();

        // copy new resources
        for (auto s : other.swarms_) {
            swarms_.push_back(s->clone());
        }

        all_members_ = other.all_members_;
        all_parents_ = other.all_parents_;
        all_parent_dists_ = other.all_parent_dists_;
        all_gens_ = other.all_gens_;
        all_rads_ = other.all_rads_;

        next_members_ = other.next_members_;
        next_parents_ = other.next_parents_;
        next_dists_ = other.next_dists_;
        next_rads_ = other.next_rads_;

        cur_members_ = other.cur_members_;
        cur_parents_ = other.cur_parents_;
        cur_dists_ = other.cur_dists_;
        cur_rads_ = other.cur_rads_;

        cur_gen_ = other.cur_gen_;
        total_num_members_ = other.total_num_members_;
        finalised_ = other.finalised_;

        return *this;

    }

    SharingSplitDacsSwarms& SharingSplitDacsSwarms::operator=(SharingSplitDacsSwarms&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        for (auto i = 0; i < swarms_.size(); i++) {
            delete swarms_[i]; swarms_[i] = nullptr;
        }
        swarms_.clear();

        // move new resources
        swarms_ = std::move(other.swarms_);
        all_members_ = std::move(other.all_members_);
        all_parents_ = std::move(other.all_parents_);
        all_parent_dists_ = std::move(other.all_parent_dists_);
        all_gens_ = std::move(other.all_gens_);
        all_rads_ = std::move(other.all_rads_);

        next_members_ = std::move(other.next_members_);
        next_parents_ = std::move(other.next_parents_);
        next_dists_ = std::move(other.next_dists_);
        next_rads_ = std::move(other.next_rads_);

        cur_members_ = std::move(other.cur_members_);
        cur_parents_ = std::move(other.cur_parents_);
        cur_dists_ = std::move(other.cur_dists_);
        cur_rads_ = std::move(other.cur_rads_);

        cur_gen_ = other.cur_gen_;
        total_num_members_ = other.total_num_members_;
        finalised_ = other.finalised_;

        return *this;

    }

    SharingSplitDacsSwarms* SharingSplitDacsSwarms::clone() const { // deep-copy clone method
        return new SharingSplitDacsSwarms(*this);
    }


    numSeqs_t SharingSplitDacsSwarms::size() const {
        return swarms_.size();
    }

    Swarm& SharingSplitDacsSwarms::get(const numSeqs_t i) {
        return *(swarms_[i]);
    }

    const Swarm& SharingSplitDacsSwarms::get(const numSeqs_t i) const {
        return *(swarms_[i]);
    }

    void SharingSplitDacsSwarms::clear() {
        std::cerr << "ERROR: Experimental SharingSplitDacsSwarms currently does not support the clear() operation." << std::endl;
    }

    void SharingSplitDacsSwarms::finalise() {

        transfer_generations();

        all_members_.finalise();
        all_parents_.finalise();
        all_parent_dists_.finalise();
        all_gens_.finalise();
        all_rads_.finalise();

        finalised_ = true;

    }

    void SharingSplitDacsSwarms::append(numSeqs_t m, numSeqs_t p, dist_t d, numSeqs_t g, dist_t r) {

        next_members_.push_back(m);
        next_parents_.push_back(p);
        next_dists_.push_back(static_cast<lenSeqs_t>(d));
        next_rads_.push_back(static_cast<lenSeqs_t>(r));

    }

    void SharingSplitDacsSwarms::swap(numSeqs_t i, numSeqs_t j) {

        std::swap(next_members_[i], next_members_[j]);
        std::swap(next_parents_[i], next_parents_[j]);
        std::swap(next_dists_[i], next_dists_[j]);
        std::swap(next_rads_[i], next_rads_[j]);

    }

    void SharingSplitDacsSwarms::sort_next_generation(const AmpliconCollection& ac) {

        std::vector<numSeqs_t> permutation(next_members_.size());
        std::iota(permutation.begin(), permutation.end(), 0);

        std::sort(permutation.begin(), permutation.end(),
                  [&](const numSeqs_t a, const numSeqs_t b) {
                      return (ac.ab(next_members_[a]) > ac.ab(next_members_[b])) ||
                             ((ac.ab(next_members_[a]) == ac.ab(next_members_[b])) &&
                              (*ac.id(next_members_[a]) < *ac.id(next_members_[b])));
                  }
        );

        std::vector<bool> done(next_members_.size());
        for (numSeqs_t i = 0; i < next_members_.size(); ++i) {

            if (!done[i]) {

                done[i] = true;
                numSeqs_t prev_j = i;
                numSeqs_t j = permutation[i];

                while (i != j) {

                    swap(prev_j, j);
                    done[j] = true;
                    prev_j = j;
                    j = permutation[j];

                }
            }

        }

        transfer_generations();

    }

    size_t SharingSplitDacsSwarms::size_in_bytes() const {

        size_t allocated_swarms = 0;
        for (auto s : swarms_) {
            allocated_swarms += s->size_in_bytes();
        }

        return sizeof(SharingSplitDacsSwarms) // SharingSplitDacsSwarms itself
            + allocated_swarms // swarms_ (content)
            + vec_pod_size_in_bytes(swarms_) - sizeof(std::vector<Swarm*>) // swarms_ (pointers, do not count size of vector object twice)
            + all_members_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) // all_members_ (do not count size of DacsCapacityArray object twice)
            + all_parents_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) // all_parents_ (do not count size of DacsCapacityArray object twice)
            + all_parent_dists_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, lenSeqs_t>) // all_parent_dists_ (do not count size of DacsCapacityArray object twice)
            + all_gens_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) // all_gens_ (do not count size of DacsCapacityArray object twice)
            + all_rads_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, lenSeqs_t>) // all_rads_ (do not count size of DacsCapacityArray object twice)
            + vec_pod_size_in_bytes(cur_members_) - sizeof(std::vector<numSeqs_t>) // cur_members_ (do not count size of vector object twice)
            + vec_pod_size_in_bytes(cur_parents_) - sizeof(std::vector<numSeqs_t>) // cur_parents_ (do not count size of vector object twice)
            + vec_pod_size_in_bytes(cur_dists_) - sizeof(std::vector<lenSeqs_t>) // cur_dists_ (do not count size of vector object twice)
            + vec_pod_size_in_bytes(cur_rads_) - sizeof(std::vector<lenSeqs_t>) // cur_rads_ (do not count size of vector object twice)
            + vec_pod_size_in_bytes(next_members_) - sizeof(std::vector<numSeqs_t>) // next_members_ (do not count size of vector object twice)
            + vec_pod_size_in_bytes(next_parents_) - sizeof(std::vector<numSeqs_t>) // next_parents_ (do not count size of vector object twice)
            + vec_pod_size_in_bytes(next_dists_) - sizeof(std::vector<lenSeqs_t>) // next_dists_ (do not count size of vector object twice)
            + vec_pod_size_in_bytes(next_rads_) - sizeof(std::vector<lenSeqs_t>); // next_rads_ (do not count size of vector object twice)

    }

    void SharingSplitDacsSwarms::show_memory(numSeqs_t pid) const {

        size_t allocated_swarms = 0;
        for (auto s : swarms_) {
            allocated_swarms += s->size_in_bytes();
        }

        std::cout << "##################################################" << std::endl;
        std::cout << "# SharingSplitDacsSwarms " << std::endl;
        std::cout << "#  - pid: " << pid << std::endl;
        std::cout << "#  - number of swarms: " << size() << std::endl;
        std::cout << "# " << std::endl;
        std::cout << "# sizeof(SharingSplitDacsSwarms): " << sizeof(SharingSplitDacsSwarms) << " bytes" << std::endl;
        std::cout << "# sizeof(Swarm*): " << sizeof(char*) << " bytes" << std::endl;
        std::cout << "# sizeof(std::vector<Swarm*>): " << sizeof(std::vector<Swarm*>) << " bytes" << std::endl;
        std::cout << "# sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>): " << sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) << " bytes" << std::endl;
        std::cout << "# sizeof(DacsCapacityArray<dacs_t, builder_t, lenSeqs_t>): " << sizeof(DacsCapacityArray<dacs_t, builder_t, lenSeqs_t>) << " bytes" << std::endl;
        std::cout << "# sizeof(std::vector<numSeqs_t>): " << sizeof(std::vector<numSeqs_t>) << " bytes" << std::endl;
        std::cout << "# sizeof(std::vector<lenSeqs_t>): " << sizeof(std::vector<lenSeqs_t>) << " bytes" << std::endl;
        std::cout << "# sizeof(numSeqs_t): " << sizeof(numSeqs_t) << " bytes" << std::endl;
        std::cout << "# sizeof(bool): " << sizeof(bool) << " bytes" << std::endl;
        std::cout << "# " << std::endl;
        std::cout << "# Total size: " << size_in_bytes() << " bytes" << std::endl;
        std::cout << "# swarms_ (content): " << allocated_swarms << " bytes" << std::endl;
        std::cout << "# swarms_ (pointers): " << vec_pod_size_in_bytes(swarms_) << " bytes" << std::endl;
        std::cout << "# all_members_: " << all_members_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# all_parents_: " << all_parents_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# all_parent_dists_: " << all_parent_dists_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# all_gens_: " << all_gens_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# all_rads_: " << all_rads_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# cur_members_: " << vec_pod_size_in_bytes(cur_members_) << " bytes" << std::endl;
        std::cout << "# cur_parents_: " << vec_pod_size_in_bytes(cur_parents_) << " bytes" << std::endl;
        std::cout << "# cur_dists_: " << vec_pod_size_in_bytes(cur_dists_) << " bytes" << std::endl;
        std::cout << "# cur_rads_: " << vec_pod_size_in_bytes(cur_rads_) << " bytes" << std::endl;
        std::cout << "# next_members_: " << vec_pod_size_in_bytes(next_members_) << " bytes" << std::endl;
        std::cout << "# next_parents_: " << vec_pod_size_in_bytes(next_parents_) << " bytes" << std::endl;
        std::cout << "# next_dists_: " << vec_pod_size_in_bytes(next_dists_) << " bytes" << std::endl;
        std::cout << "# next_rads_: " << vec_pod_size_in_bytes(next_rads_) << " bytes" << std::endl;
        std::cout << "##################################################" << std::endl;

    }

    void SharingSplitDacsSwarms::transfer_generations() {

        for (auto i = 0; i < cur_members_.size(); i++) {

            all_members_.push_back(cur_members_[i]);
            all_parents_.push_back(cur_parents_[i]);
            all_parent_dists_.push_back(cur_dists_[i]);
            all_gens_.push_back(cur_gen_);
            all_rads_.push_back(cur_rads_[i]);

        }
        total_num_members_ += cur_members_.size();

        cur_members_.swap(next_members_);
        cur_parents_.swap(next_parents_);
        cur_dists_.swap(next_dists_);
        cur_rads_.swap(next_rads_);
        cur_gen_++;

        std::vector<numSeqs_t>().swap(next_members_);
        std::vector<numSeqs_t>().swap(next_parents_);
        std::vector<lenSeqs_t>().swap(next_dists_);
        std::vector<lenSeqs_t>().swap(next_rads_);

    }


    numSeqs_t SharingSplitDacsSwarms::member(numSeqs_t i) const {

        if (finalised_) return all_members_.get(i);
        auto pos = i - total_num_members_;
        return (pos < cur_members_.size()) ? cur_members_[pos] : next_members_[pos - cur_members_.size()];

    }

    numSeqs_t SharingSplitDacsSwarms::parent(numSeqs_t i) const {

        if (finalised_) return all_parents_.get(i);
        auto pos = i - total_num_members_;
        return (pos < cur_parents_.size()) ? cur_parents_[pos] : next_parents_[pos - cur_parents_.size()];

    }

    numSeqs_t SharingSplitDacsSwarms::gen(numSeqs_t i) const {

        if (finalised_) return all_gens_.get(i);
        auto pos = i - total_num_members_;
        return (pos < cur_members_.size()) ? cur_gen_ : (cur_gen_ + 1);

    }

    dist_t SharingSplitDacsSwarms::rad(numSeqs_t i) const {

        if (finalised_) return all_rads_.get(i);
        auto pos = i - total_num_members_;
        return (pos < cur_rads_.size()) ? cur_rads_[pos] : next_rads_[pos - cur_rads_.size()];

    }

    dist_t SharingSplitDacsSwarms::parent_dist(numSeqs_t i) const {

        if (finalised_) return all_parent_dists_.get(i);
        auto pos = i - total_num_members_;
        return (pos < cur_dists_.size()) ? cur_dists_[pos] : next_dists_[pos - cur_dists_.size()];

    }


    /* === SharingSplitDacsSwarms::ForwardingDacsSwarm === */

    SharingSplitDacsSwarms::ForwardingDacsSwarm::ForwardingDacsSwarm(SharingSplitDacsSwarms& as, numSeqs_t ab) {

        swarms_ = &as;
        start_ = as.total_num_members_;

        size_ = 1;
        mass_ = ab;
        num_different_ = 1;
        num_singletons_ = (ab == 1);

        max_rad_ = 0;

        graftings_ = nullptr;
        attached_ = false;

    }

    SharingSplitDacsSwarms::ForwardingDacsSwarm::~ForwardingDacsSwarm() {
        delete graftings_;
    }

    SharingSplitDacsSwarms::ForwardingDacsSwarm::ForwardingDacsSwarm(const ForwardingDacsSwarm& other) : // copy constructor
            swarms_(other.swarms_), start_(other.start_), size_(other.size_), mass_(other.mass_), num_different_(other.num_different_),
            num_singletons_(other.num_singletons_), max_rad_(other.max_rad_),
            graftings_(other.graftings_ ? other.graftings_->clone() : nullptr), attached_(other.attached_) {

        // nothing else to do

    }

    SharingSplitDacsSwarms::ForwardingDacsSwarm::ForwardingDacsSwarm(ForwardingDacsSwarm&& other) noexcept : // move constructor
            swarms_(other.swarms_), start_(other.start_),
            size_(other.size_), mass_(other.mass_), num_different_(other.num_different_), num_singletons_(other.num_singletons_),
            max_rad_(other.max_rad_), graftings_(other.graftings_), attached_(other.attached_) {

        other.graftings_ = nullptr;

    }

    SharingSplitDacsSwarms::ForwardingDacsSwarm& SharingSplitDacsSwarms::ForwardingDacsSwarm::operator=(const ForwardingDacsSwarm& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete graftings_;

        // copy new resources
        swarms_ = other.swarms_;
        start_ = other.start_;

        size_ = other.size_;
        mass_ = other.mass_;
        num_different_ = other.num_different_;
        num_singletons_ = other.num_singletons_;

        max_rad_ = other.max_rad_;

        graftings_ = other.graftings_ ? other.graftings_->clone() : nullptr;
        attached_ = other.attached_;

        return *this;

    }

    SharingSplitDacsSwarms::ForwardingDacsSwarm& SharingSplitDacsSwarms::ForwardingDacsSwarm::operator=(ForwardingDacsSwarm&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete graftings_;

        // copy new resources
        swarms_ = other.swarms_;
        start_ = other.start_;

        size_ = other.size_;
        mass_ = other.mass_;
        num_different_ = other.num_different_;
        num_singletons_ = other.num_singletons_;

        max_rad_ = other.max_rad_;

        graftings_ = other.graftings_; other.graftings_ = nullptr;
        attached_ = other.attached_;

        return *this;

    }

    SharingSplitDacsSwarms::ForwardingDacsSwarm* SharingSplitDacsSwarms::ForwardingDacsSwarm::clone() const {
        return new ForwardingDacsSwarm(*this);
    }


    numSeqs_t SharingSplitDacsSwarms::ForwardingDacsSwarm::seed() const {
        return swarms_->member(start_);
    }
    SwarmEntry* SharingSplitDacsSwarms::ForwardingDacsSwarm::seed_entry() const {
        return new ForwardingSwarmEntry<ForwardingDacsSwarm>(this, 0);
    }
    SwarmEntry* SharingSplitDacsSwarms::ForwardingDacsSwarm::seed_entry(SwarmEntry* entry) const {

        dynamic_cast<ForwardingSwarmEntry<ForwardingDacsSwarm>*>(entry)->reassign(this, 0);
        return entry;

    }

    numSeqs_t SharingSplitDacsSwarms::ForwardingDacsSwarm::get(numSeqs_t i) {
        return swarms_->member(start_ + i);
    }
    SwarmEntry* SharingSplitDacsSwarms::ForwardingDacsSwarm::get_entry(numSeqs_t i) {
        return new ForwardingSwarmEntry<ForwardingDacsSwarm>(this, i);
    }
    SwarmEntry* SharingSplitDacsSwarms::ForwardingDacsSwarm::get_entry(numSeqs_t i, SwarmEntry* entry) const {

        dynamic_cast<ForwardingSwarmEntry<ForwardingDacsSwarm>*>(entry)->reassign(this, i);
        return entry;

    }

    numSeqs_t SharingSplitDacsSwarms::ForwardingDacsSwarm::member(numSeqs_t i) const {
        return swarms_->member(start_ + i);
    }

    numSeqs_t SharingSplitDacsSwarms::ForwardingDacsSwarm::parent(numSeqs_t i) const {
        return swarms_->parent(start_ + i);
    }

    numSeqs_t SharingSplitDacsSwarms::ForwardingDacsSwarm::gen(numSeqs_t i) const {
        return swarms_->gen(start_ + i);
    }

    dist_t SharingSplitDacsSwarms::ForwardingDacsSwarm::rad(numSeqs_t i) const {
        return swarms_->rad(start_ + i);
    }

    dist_t SharingSplitDacsSwarms::ForwardingDacsSwarm::parent_dist(numSeqs_t i) const {
        return swarms_->parent_dist(start_ + i);
    }

    numSeqs_t SharingSplitDacsSwarms::ForwardingDacsSwarm::size() const {
        return size_;
    }

    numSeqs_t SharingSplitDacsSwarms::ForwardingDacsSwarm::num_different() const {
        return num_different_;
    }

    numSeqs_t SharingSplitDacsSwarms::ForwardingDacsSwarm::num_singletons() const {
        return num_singletons_;
    }

    numSeqs_t SharingSplitDacsSwarms::ForwardingDacsSwarm::mass() const {
        return mass_;
    }

    void SharingSplitDacsSwarms::ForwardingDacsSwarm::sort_next_generation(const AmpliconCollection& ac, numSeqs_t p) {
        swarms_->sort_next_generation(ac);
    }

    void SharingSplitDacsSwarms::ForwardingDacsSwarm::append(numSeqs_t m, numSeqs_t g, dist_t r, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) {

        swarms_->append(m, p, dist, g, r);

        mass_ += ab;
        num_different_ += new_seq;
        num_singletons_ += (ab == 1);

        max_rad_ = std::max(max_rad_, r);

        size_++;

    }

    void SharingSplitDacsSwarms::ForwardingDacsSwarm::append(numSeqs_t m, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) {

        numSeqs_t cur_pos = start_ + p - swarms_->total_num_members_;
        swarms_->append(m, swarms_->cur_members_[cur_pos], dist, swarms_->cur_gen_ + 1, swarms_->cur_rads_[cur_pos] + dist);

        mass_ += ab;
        num_different_ += new_seq;
        num_singletons_ += (ab == 1);

        max_rad_ = std::max(max_rad_, swarms_->cur_rads_[cur_pos] + dist);

        size_++;

    }

    void SharingSplitDacsSwarms::ForwardingDacsSwarm::attach(const numSeqs_t p, const numSeqs_t c, Swarm* cs, const dist_t d) {

        if (graftings_ == nullptr) {
            graftings_ = new SimpleGraftingInfo();
        }

        graftings_->add_link(p, c, cs, d);
        cs->mark_as_attached();

    }

    void SharingSplitDacsSwarms::ForwardingDacsSwarm::mark_as_attached() {
        attached_ = true;
    }

    bool SharingSplitDacsSwarms::ForwardingDacsSwarm::is_attached() const {
        return attached_;
    }

    numSeqs_t SharingSplitDacsSwarms::ForwardingDacsSwarm::total_size() const {

        numSeqs_t sum = size();
        if (graftings_ != nullptr) {

            for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                sum += graftings_->get_child_swarm(i)->size();
            }

        }

        return sum;

    }

    numSeqs_t SharingSplitDacsSwarms::ForwardingDacsSwarm::total_num_different() const {

        numSeqs_t sum = num_different();
        if (graftings_ != nullptr) {

            for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                sum += graftings_->get_child_swarm(i)->num_different();
            }

        }

        return sum;

    }

    numSeqs_t SharingSplitDacsSwarms::ForwardingDacsSwarm::total_num_singletons() const {

        numSeqs_t sum = num_singletons();
        if (graftings_ != nullptr) {

            for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                sum += graftings_->get_child_swarm(i)->num_singletons();
            }

        }

        return sum;

    }

    numSeqs_t SharingSplitDacsSwarms::ForwardingDacsSwarm::total_mass() const {

        numSeqs_t sum = mass();
        if (graftings_ != nullptr) {

            for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                sum += graftings_->get_child_swarm(i)->mass();
            }

        }

        return sum;

    }

    dist_t SharingSplitDacsSwarms::ForwardingDacsSwarm::max_rad() const {
        return max_rad_;
    }

    lenSeqs_t SharingSplitDacsSwarms::ForwardingDacsSwarm::max_gen() const {
        return (size_ > 0) ? swarms_->gen(start_ + size_ - 1) : 0;
    }

    GraftingInfo* SharingSplitDacsSwarms::ForwardingDacsSwarm::get_grafting_info() const {
        return graftings_;
    }

    size_t SharingSplitDacsSwarms::ForwardingDacsSwarm::size_in_bytes() const {
        return sizeof(ForwardingDacsSwarm) // ForwardingDacsSwarm itself (no ownership over swarms_)
            + ((graftings_ != nullptr) ? graftings_->size_in_bytes() : 0); // graftings_
    }



    /* === SharingSplitFullDacsSwarms === */

    SharingSplitFullDacsSwarms::SharingSplitFullDacsSwarms(const numSeqs_t pool_size, const ChunkLengths& config)
            : all_members_(config.ids, pool_size, config.level_factor, config.extend_factor, config.opt),
              all_parents_(config.ids, pool_size, config.level_factor, config.extend_factor, config.opt),
              all_parent_dists_(config.dists, pool_size, config.level_factor, config.extend_factor, config.opt),
              all_gens_(config.gens, pool_size, config.level_factor, config.extend_factor, config.opt),
              all_rads_(config.rads, pool_size, config.level_factor, config.extend_factor, config.opt),
              swarms_start_(config.starts, config.get_capa(pool_size), config.level_factor, config.extend_factor, config.opt),
              swarms_size_(config.sizes, config.get_capa(pool_size), config.level_factor, config.extend_factor, config.opt),
              swarms_mass_(config.masses, config.get_capa(pool_size), config.level_factor, config.extend_factor, config.opt),
              swarms_num_different_(config.differents, config.get_capa(pool_size), config.level_factor, config.extend_factor, config.opt),
              swarms_num_singletons_(config.singletons, config.get_capa(pool_size), config.level_factor, config.extend_factor, config.opt),
              swarms_max_rad_(config.max_rads, config.get_capa(pool_size), config.level_factor, config.extend_factor, config.opt) {

        cur_gen_ = 0;
        total_num_members_ = 0;
        finalised_ = false;

        cur_swarm_start_ = 0;
        cur_swarm_size_ = 0;
        cur_swarm_mass_ = 0;
        cur_swarm_num_different_ = 0;
        cur_swarm_num_singletons_ = 0;
        cur_swarm_max_rad_ = 0;

    }

    Swarm& SharingSplitFullDacsSwarms::initialise_cluster(numSeqs_t s, numSeqs_t ab) {

        transfer_generations();
        close_swarm(s, ab);

        swarms_.push_back(new ForwardingDacsSwarm(this, swarms_.size()));
        cur_swarm_start_ = total_num_members_;
        cur_swarm_size_ = 1;
        cur_swarm_mass_ = ab;
        cur_swarm_num_different_ = 1;
        cur_swarm_num_singletons_ = (ab == 1);
        cur_swarm_max_rad_ = 0;

        cur_members_.push_back(s);
        cur_parents_.push_back(s);
        cur_dists_.push_back(0);
        cur_gen_ = 0;
        cur_rads_.push_back(0);

        return *(swarms_.back());

    }

    SharingSplitFullDacsSwarms::~SharingSplitFullDacsSwarms() {

        for (auto i = 0; i < swarms_.size(); i++) {
            delete swarms_[i]; swarms_[i] = nullptr;
        }
        for (auto i = 0; i < swarms_graftings_.size(); i++) {
            delete swarms_graftings_[i]; swarms_graftings_[i] = nullptr;
        }

    }

    SharingSplitFullDacsSwarms::SharingSplitFullDacsSwarms(const SharingSplitFullDacsSwarms& other) :
            swarms_start_(other.swarms_start_), swarms_size_(other.swarms_size_), swarms_mass_(other.swarms_mass_),
            swarms_num_different_(other.swarms_num_different_), swarms_num_singletons_(other.swarms_num_singletons_),
            swarms_max_rad_(other.swarms_max_rad_), swarms_attached_(other.swarms_attached_),
            all_members_(other.all_members_), all_parents_(other.all_parents_),
            all_parent_dists_(other.all_parent_dists_), all_gens_(other.all_gens_), all_rads_(other.all_rads_),
            next_members_(other.next_members_), next_parents_(other.next_parents_),
            next_dists_(other.next_dists_), next_rads_(other.next_rads_),
            cur_members_(other.cur_members_), cur_parents_(other.cur_parents_),
            cur_dists_(other.cur_dists_), cur_rads_(other.cur_rads_) { // copy constructor

        cur_swarm_start_ = other.cur_swarm_start_;
        cur_swarm_size_ = other.cur_swarm_size_;
        cur_swarm_mass_ = other.cur_swarm_mass_;
        cur_swarm_num_different_ = other.cur_swarm_num_different_;
        cur_swarm_num_singletons_ = other.cur_swarm_num_singletons_;
        cur_swarm_max_rad_ = other.cur_swarm_max_rad_;

        cur_gen_ = other.cur_gen_;
        total_num_members_ = other.total_num_members_;
        finalised_ = other.finalised_;

        swarms_.reserve(other.swarms_.size());
        for (auto s : other.swarms_) {
            swarms_.push_back(s->clone());
        }
        swarms_graftings_.reserve(other.swarms_graftings_.size());
        for (auto g : other.swarms_graftings_) {
            swarms_graftings_.push_back(g ? g->clone() : nullptr);
        }

    }

    SharingSplitFullDacsSwarms::SharingSplitFullDacsSwarms(SharingSplitFullDacsSwarms&& other) noexcept :
            swarms_start_(std::move(other.swarms_start_)), swarms_size_(std::move(other.swarms_size_)), swarms_mass_(std::move(other.swarms_mass_)),
            swarms_num_different_(std::move(other.swarms_num_different_)), swarms_num_singletons_(std::move(other.swarms_num_singletons_)),
            swarms_max_rad_(std::move(other.swarms_max_rad_)), swarms_attached_(std::move(other.swarms_attached_)),
            all_members_(std::move(other.all_members_)), all_parents_(std::move(other.all_parents_)),
            all_parent_dists_(std::move(other.all_parent_dists_)), all_gens_(std::move(other.all_gens_)), all_rads_(std::move(other.all_rads_)),
            next_members_(std::move(other.next_members_)), next_parents_(std::move(other.next_parents_)),
            next_dists_(std::move(other.next_dists_)), next_rads_(std::move(other.next_rads_)),
            cur_members_(std::move(other.cur_members_)), cur_parents_(std::move(other.cur_parents_)),
            cur_dists_(std::move(other.cur_dists_)), cur_rads_(std::move(other.cur_rads_)) { // move constructor

        cur_swarm_start_ = other.cur_swarm_start_;
        cur_swarm_size_ = other.cur_swarm_size_;
        cur_swarm_mass_ = other.cur_swarm_mass_;
        cur_swarm_num_different_ = other.cur_swarm_num_different_;
        cur_swarm_num_singletons_ = other.cur_swarm_num_singletons_;
        cur_swarm_max_rad_ = other.cur_swarm_max_rad_;

        cur_gen_ = other.cur_gen_;
        total_num_members_ = other.total_num_members_;
        finalised_ = other.finalised_;

        swarms_.reserve(other.swarms_.size());
        for (auto i = 0; i < other.swarms_.size(); i++) {
            swarms_.push_back(other.swarms_[i]); other.swarms_[i] = nullptr;
        }
        swarms_graftings_.reserve(other.swarms_graftings_.size());
        for (auto i = 0; i < other.swarms_graftings_.size(); i++) {
            swarms_graftings_.push_back(other.swarms_graftings_[i]); other.swarms_graftings_[i] = nullptr;
        }

    }

    SharingSplitFullDacsSwarms& SharingSplitFullDacsSwarms::operator=(const SharingSplitFullDacsSwarms& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        for (auto i = 0; i < swarms_.size(); i++) {
            delete swarms_[i]; swarms_[i] = nullptr;
        }
        swarms_.clear();
        for (auto i = 0; i < swarms_graftings_.size(); i++) {
            delete swarms_graftings_[i]; swarms_graftings_[i] = nullptr;
        }
        swarms_graftings_.clear();

        // copy new resources
        swarms_.reserve(other.swarms_.size());
        for (auto s : other.swarms_) {
            swarms_.push_back(s->clone());
        }
        swarms_graftings_.reserve(other.swarms_graftings_.size());
        for (auto g : other.swarms_graftings_) {
            swarms_graftings_.push_back(g ? g->clone() : nullptr);
        }

        swarms_start_ = other.swarms_start_;
        swarms_size_ = other.swarms_size_;
        swarms_mass_ = other.swarms_mass_;
        swarms_num_different_ = other.swarms_num_different_;
        swarms_num_singletons_ = other.swarms_num_singletons_;
        swarms_max_rad_ = other.swarms_max_rad_;
        swarms_attached_ = other.swarms_attached_;
        all_members_ = other.all_members_;
        all_parents_ = other.all_parents_;
        all_parent_dists_ = other.all_parent_dists_;
        all_gens_ = other.all_gens_;
        all_rads_ = other.all_rads_;
        next_members_ = other.next_members_;
        next_parents_ = other.next_parents_;
        next_dists_ = other.next_dists_;
        next_rads_ = other.next_rads_;
        cur_members_ = other.cur_members_;
        cur_parents_ = other.cur_parents_;
        cur_dists_ = other.cur_dists_;
        cur_rads_ = other.cur_rads_;

        cur_swarm_start_ = other.cur_swarm_start_;
        cur_swarm_size_ = other.cur_swarm_size_;
        cur_swarm_mass_ = other.cur_swarm_mass_;
        cur_swarm_num_different_ = other.cur_swarm_num_different_;
        cur_swarm_num_singletons_ = other.cur_swarm_num_singletons_;
        cur_swarm_max_rad_ = other.cur_swarm_max_rad_;

        cur_gen_ = other.cur_gen_;
        total_num_members_ = other.total_num_members_;
        finalised_ = other.finalised_;

        return *this;

    }

    SharingSplitFullDacsSwarms& SharingSplitFullDacsSwarms::operator=(SharingSplitFullDacsSwarms&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        for (auto i = 0; i < swarms_.size(); i++) {
            delete swarms_[i]; swarms_[i] = nullptr;
        }
        swarms_.clear();
        for (auto i = 0; i < swarms_graftings_.size(); i++) {
            delete swarms_graftings_[i]; swarms_graftings_[i] = nullptr;
        }
        swarms_graftings_.clear();

        // move new resources
        swarms_.reserve(other.swarms_.size());
        for (auto i = 0; i < other.swarms_.size(); i++) {
            swarms_.push_back(other.swarms_[i]); other.swarms_[i] = nullptr;
        }
        swarms_graftings_.reserve(other.swarms_graftings_.size());
        for (auto i = 0; i < other.swarms_graftings_.size(); i++) {
            swarms_graftings_.push_back(other.swarms_graftings_[i]); other.swarms_graftings_[i] = nullptr;
        }

        swarms_start_ = std::move(other.swarms_start_);
        swarms_size_ = std::move(other.swarms_size_);
        swarms_mass_ = std::move(other.swarms_mass_);
        swarms_num_different_ = std::move(other.swarms_num_different_);
        swarms_num_singletons_ = std::move(other.swarms_num_singletons_);
        swarms_max_rad_ = std::move(other.swarms_max_rad_);
        swarms_attached_ = std::move(other.swarms_attached_);
        all_members_ = std::move(other.all_members_);
        all_parents_ = std::move(other.all_parents_);
        all_parent_dists_ = std::move(other.all_parent_dists_);
        all_gens_ = std::move(other.all_gens_);
        all_rads_ = std::move(other.all_rads_);
        next_members_ = std::move(other.next_members_);
        next_parents_ = std::move(other.next_parents_);
        next_dists_ = std::move(other.next_dists_);
        next_rads_ = std::move(other.next_rads_);
        cur_members_ = std::move(other.cur_members_);
        cur_parents_ = std::move(other.cur_parents_);
        cur_dists_ = std::move(other.cur_dists_);
        cur_rads_ = std::move(other.cur_rads_);

        cur_swarm_start_ = other.cur_swarm_start_;
        cur_swarm_size_ = other.cur_swarm_size_;
        cur_swarm_mass_ = other.cur_swarm_mass_;
        cur_swarm_num_different_ = other.cur_swarm_num_different_;
        cur_swarm_num_singletons_ = other.cur_swarm_num_singletons_;
        cur_swarm_max_rad_ = other.cur_swarm_max_rad_;

        cur_gen_ = other.cur_gen_;
        total_num_members_ = other.total_num_members_;
        finalised_ = other.finalised_;

        return *this;

    }

    SharingSplitFullDacsSwarms* SharingSplitFullDacsSwarms::clone() const { // deep-copy clone method
        return new SharingSplitFullDacsSwarms(*this);
    }


    numSeqs_t SharingSplitFullDacsSwarms::size() const {
        return swarms_.size();
    }

    Swarm& SharingSplitFullDacsSwarms::get(const numSeqs_t i) {
        return *(swarms_[i]);
    }

    const Swarm& SharingSplitFullDacsSwarms::get(const numSeqs_t i) const {
        return *(swarms_[i]);
    }

    void SharingSplitFullDacsSwarms::clear() {
        std::cerr << "ERROR: Experimental SharingSplitFullDacsSwarms currently does not support the clear() operation." << std::endl;
    }

    void SharingSplitFullDacsSwarms::finalise() {

        transfer_generations();
        close_swarm(0, 0);

        all_members_.finalise();
        all_parents_.finalise();
        all_parent_dists_.finalise();
        all_gens_.finalise();
        all_rads_.finalise();

        swarms_start_.finalise();
        swarms_size_.finalise();
        swarms_mass_.finalise();
        swarms_num_different_.finalise();
        swarms_num_singletons_.finalise();
        swarms_max_rad_.finalise();

        finalised_ = true;

    }

    void SharingSplitFullDacsSwarms::append(numSeqs_t m, numSeqs_t p, dist_t d, numSeqs_t g, dist_t r, numSeqs_t ab, bool new_seq) {

        next_members_.push_back(m);
        next_parents_.push_back(p);
        next_dists_.push_back(static_cast<lenSeqs_t>(d));
        next_rads_.push_back(static_cast<lenSeqs_t>(r));

        cur_swarm_mass_ += ab;
        cur_swarm_num_different_ += new_seq;
        cur_swarm_num_singletons_ += (ab == 1);
        cur_swarm_max_rad_ = std::max(cur_swarm_max_rad_, r);
        cur_swarm_size_++;

    }

    void SharingSplitFullDacsSwarms::swap(numSeqs_t i, numSeqs_t j) {

        std::swap(next_members_[i], next_members_[j]);
        std::swap(next_parents_[i], next_parents_[j]);
        std::swap(next_dists_[i], next_dists_[j]);
        std::swap(next_rads_[i], next_rads_[j]);

    }

    void SharingSplitFullDacsSwarms::sort_next_generation(const AmpliconCollection& ac) {

        std::vector<numSeqs_t> permutation(next_members_.size());
        std::iota(permutation.begin(), permutation.end(), 0);

        std::sort(permutation.begin(), permutation.end(),
                  [&](const numSeqs_t a, const numSeqs_t b) {
                      return (ac.ab(next_members_[a]) > ac.ab(next_members_[b])) ||
                             ((ac.ab(next_members_[a]) == ac.ab(next_members_[b])) &&
                              (*ac.id(next_members_[a]) < *ac.id(next_members_[b])));
                  }
        );

        std::vector<bool> done(next_members_.size());
        for (numSeqs_t i = 0; i < next_members_.size(); ++i) {

            if (!done[i]) {

                done[i] = true;
                numSeqs_t prev_j = i;
                numSeqs_t j = permutation[i];

                while (i != j) {

                    swap(prev_j, j);
                    done[j] = true;
                    prev_j = j;
                    j = permutation[j];

                }
            }

        }

        transfer_generations();

    }

    size_t SharingSplitFullDacsSwarms::size_in_bytes() const {

        size_t allocated_swarms = 0;
        for (auto s : swarms_) {
            allocated_swarms += s->size_in_bytes();
        }
        size_t allocated_graftings = 0;
        for (auto g : swarms_graftings_) {
            allocated_graftings += (g != nullptr) ? g->size_in_bytes() : 0;
        }

        return sizeof(SharingSplitFullDacsSwarms) // SharingSplitFullDacsSwarms itself
               + allocated_swarms // swarms_ (content)
               + vec_pod_size_in_bytes(swarms_) - sizeof(std::vector<Swarm*>) // swarms_ (pointers, do not count size of vector object twice)
               + all_members_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) // all_members_ (do not count size of DacsCapacityArray object twice)
               + all_parents_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) // all_parents_ (do not count size of DacsCapacityArray object twice)
               + all_parent_dists_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, lenSeqs_t>) // all_parent_dists_ (do not count size of DacsCapacityArray object twice)
               + all_gens_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) // all_gens_ (do not count size of DacsCapacityArray object twice)
               + all_rads_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, lenSeqs_t>) // all_rads_ (do not count size of DacsCapacityArray object twice)
               + vec_pod_size_in_bytes(cur_members_) - sizeof(std::vector<numSeqs_t>) // cur_members_ (do not count size of vector object twice)
               + vec_pod_size_in_bytes(cur_parents_) - sizeof(std::vector<numSeqs_t>) // cur_parents_ (do not count size of vector object twice)
               + vec_pod_size_in_bytes(cur_dists_) - sizeof(std::vector<lenSeqs_t>) // cur_dists_ (do not count size of vector object twice)
               + vec_pod_size_in_bytes(cur_rads_) - sizeof(std::vector<lenSeqs_t>) // cur_rads_ (do not count size of vector object twice)
               + vec_pod_size_in_bytes(next_members_) - sizeof(std::vector<numSeqs_t>) // next_members_ (do not count size of vector object twice)
               + vec_pod_size_in_bytes(next_parents_) - sizeof(std::vector<numSeqs_t>) // next_parents_ (do not count size of vector object twice)
               + vec_pod_size_in_bytes(next_dists_) - sizeof(std::vector<lenSeqs_t>) // next_dists_ (do not count size of vector object twice)
               + vec_pod_size_in_bytes(next_rads_) - sizeof(std::vector<lenSeqs_t>) // next_rads_ (do not count size of vector object twice)
               + swarms_start_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) // swarms_start_
               + swarms_size_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) // swarms_size_
               + swarms_mass_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) // swarms_mass_
               + swarms_num_different_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) // swarms_num_different_
               + swarms_num_singletons_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) // swarms_num_singletons_
               + swarms_max_rad_.size_in_bytes() - sizeof(DacsCapacityArray<dacs_t, builder_t, lenSeqs_t>) // swarms_max_rad_
               + allocated_graftings // swarms_graftings_ (content)
               + vec_pod_size_in_bytes(swarms_graftings_) - sizeof(std::vector<GraftingInfo*>) // swarms_graftings_ (pointers)
               + vec_bool_size_in_bytes(swarms_attached_) - sizeof(std::vector<bool>); // swarms_attached_

    }

    void SharingSplitFullDacsSwarms::show_memory(numSeqs_t pid) const {

        size_t allocated_swarms = 0;
        for (auto s : swarms_) {
            allocated_swarms += s->size_in_bytes();
        }
        size_t allocated_graftings = 0;
        for (auto g : swarms_graftings_) {
            allocated_graftings += (g != nullptr) ? g->size_in_bytes() : 0;
        }

        std::cout << "##################################################" << std::endl;
        std::cout << "# SharingSplitFullDacsSwarms " << std::endl;
        std::cout << "#  - pid: " << pid << std::endl;
        std::cout << "#  - number of swarms: " << size() << std::endl;
        std::cout << "# " << std::endl;
        std::cout << "# sizeof(SharingSplitFullDacsSwarms): " << sizeof(SharingSplitFullDacsSwarms) << " bytes" << std::endl;
        std::cout << "# sizeof(Swarm*): " << sizeof(char*) << " bytes" << std::endl;
        std::cout << "# sizeof(std::vector<Swarm*>): " << sizeof(std::vector<Swarm*>) << " bytes" << std::endl;
        std::cout << "# sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>): " << sizeof(DacsCapacityArray<dacs_t, builder_t, numSeqs_t>) << " bytes" << std::endl;
        std::cout << "# sizeof(DacsCapacityArray<dacs_t, builder_t, lenSeqs_t>): " << sizeof(DacsCapacityArray<dacs_t, builder_t, lenSeqs_t>) << " bytes" << std::endl;
        std::cout << "# sizeof(std::vector<numSeqs_t>): " << sizeof(std::vector<numSeqs_t>) << " bytes" << std::endl;
        std::cout << "# sizeof(std::vector<lenSeqs_t>): " << sizeof(std::vector<lenSeqs_t>) << " bytes" << std::endl;
        std::cout << "# sizeof(numSeqs_t): " << sizeof(numSeqs_t) << " bytes" << std::endl;
        std::cout << "# sizeof(bool): " << sizeof(bool) << " bytes" << std::endl;
        std::cout << "# " << std::endl;
        std::cout << "# Total size: " << size_in_bytes() << " bytes" << std::endl;
        std::cout << "# swarms_ (content): " << allocated_swarms << " bytes" << std::endl;
        std::cout << "# swarms_ (pointers): " << vec_pod_size_in_bytes(swarms_) << " bytes" << std::endl;
        std::cout << "# all_members_: " << all_members_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# all_parents_: " << all_parents_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# all_parent_dists_: " << all_parent_dists_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# all_gens_: " << all_gens_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# all_rads_: " << all_rads_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# cur_members_: " << vec_pod_size_in_bytes(cur_members_) << " bytes" << std::endl;
        std::cout << "# cur_parents_: " << vec_pod_size_in_bytes(cur_parents_) << " bytes" << std::endl;
        std::cout << "# cur_dists_: " << vec_pod_size_in_bytes(cur_dists_) << " bytes" << std::endl;
        std::cout << "# cur_rads_: " << vec_pod_size_in_bytes(cur_rads_) << " bytes" << std::endl;
        std::cout << "# next_members_: " << vec_pod_size_in_bytes(next_members_) << " bytes" << std::endl;
        std::cout << "# next_parents_: " << vec_pod_size_in_bytes(next_parents_) << " bytes" << std::endl;
        std::cout << "# next_dists_: " << vec_pod_size_in_bytes(next_dists_) << " bytes" << std::endl;
        std::cout << "# next_rads_: " << vec_pod_size_in_bytes(next_rads_) << " bytes" << std::endl;
        std::cout << "# swarms_start_: " << swarms_start_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# swarms_size_: " << swarms_size_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# swarms_mass_: " << swarms_mass_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# swarms_num_different_: " << swarms_num_different_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# swarms_num_singletons_: " << swarms_num_singletons_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# swarms_max_rad_: " << swarms_max_rad_.size_in_bytes() << " bytes" << std::endl;
        std::cout << "# swarms_graftings_ (content): " << allocated_graftings << " bytes" << std::endl;
        std::cout << "# swarms_graftings_ (pointers): " << vec_pod_size_in_bytes(swarms_graftings_) << " bytes" << std::endl;
        std::cout << "# swarms_attached_: " << vec_bool_size_in_bytes(swarms_attached_) << " bytes" << std::endl;
        std::cout << "##################################################" << std::endl;

    }

    void SharingSplitFullDacsSwarms::transfer_generations() {

        for (auto i = 0; i < cur_members_.size(); i++) {

            all_members_.push_back(cur_members_[i]);
            all_parents_.push_back(cur_parents_[i]);
            all_parent_dists_.push_back(cur_dists_[i]);
            all_gens_.push_back(cur_gen_);
            all_rads_.push_back(cur_rads_[i]);

        }
        total_num_members_ += cur_members_.size();

        cur_members_.swap(next_members_);
        cur_parents_.swap(next_parents_);
        cur_dists_.swap(next_dists_);
        cur_rads_.swap(next_rads_);
        cur_gen_++;

        std::vector<numSeqs_t>().swap(next_members_);
        std::vector<numSeqs_t>().swap(next_parents_);
        std::vector<lenSeqs_t>().swap(next_dists_);
        std::vector<lenSeqs_t>().swap(next_rads_);

    }

    void SharingSplitFullDacsSwarms::close_swarm(numSeqs_t s, numSeqs_t ab) {

        if (cur_swarm_mass_ > 0) {

            swarms_start_.push_back(cur_swarm_start_);
            swarms_size_.push_back(cur_swarm_size_);
            swarms_mass_.push_back(cur_swarm_mass_);
            swarms_num_different_.push_back(cur_swarm_num_different_);
            swarms_num_singletons_.push_back(cur_swarm_num_singletons_);
            swarms_max_rad_.push_back(cur_swarm_max_rad_);
            swarms_graftings_.push_back(nullptr);
            swarms_attached_.push_back(false);

        }

        cur_swarm_start_ = total_num_members_;
        cur_swarm_size_ = 1;
        cur_swarm_mass_ = ab;
        cur_swarm_num_different_ = 1;
        cur_swarm_num_singletons_ = (ab == 1);
        cur_swarm_max_rad_ = 0;

    }


    numSeqs_t SharingSplitFullDacsSwarms::member(numSeqs_t sid, numSeqs_t i) const {

        if (finalised_) return all_members_.get(swarms_start_.get(sid) + i);
        auto pos = cur_swarm_start_ + i - total_num_members_;
        return (pos < cur_members_.size()) ? cur_members_[pos] : next_members_[pos - cur_members_.size()];

    }

    numSeqs_t SharingSplitFullDacsSwarms::parent(numSeqs_t sid, numSeqs_t i) const {

        if (finalised_) return all_parents_.get(swarms_start_.get(sid) + i);
        auto pos = cur_swarm_start_ + i - total_num_members_;
        return (pos < cur_parents_.size()) ? cur_parents_[pos] : next_parents_[pos - cur_parents_.size()];

    }

    numSeqs_t SharingSplitFullDacsSwarms::gen(numSeqs_t sid, numSeqs_t i) const {

        if (finalised_) return all_gens_.get(swarms_start_.get(sid) + i);
        auto pos = cur_swarm_start_ + i - total_num_members_;
        return (pos < cur_members_.size()) ? cur_gen_ : (cur_gen_ + 1);

    }

    dist_t SharingSplitFullDacsSwarms::rad(numSeqs_t sid, numSeqs_t i) const {

        if (finalised_) return all_rads_.get(swarms_start_.get(sid) + i);
        auto pos = cur_swarm_start_ + i - total_num_members_;
        return (pos < cur_rads_.size()) ? cur_rads_[pos] : next_rads_[pos - cur_rads_.size()];

    }

    dist_t SharingSplitFullDacsSwarms::parent_dist(numSeqs_t sid, numSeqs_t i) const {

        if (finalised_) return all_parent_dists_.get(swarms_start_.get(sid) + i);
        auto pos = cur_swarm_start_ + i - total_num_members_;
        return (pos < cur_dists_.size()) ? cur_dists_[pos] : next_dists_[pos - cur_dists_.size()];

    }

    numSeqs_t SharingSplitFullDacsSwarms::swarm_size(numSeqs_t sid) const {
        return (finalised_) ? swarms_size_.get(sid) : cur_swarm_size_;
    }

    numSeqs_t SharingSplitFullDacsSwarms::swarm_num_different(numSeqs_t sid) const {
        return (finalised_) ? swarms_num_different_.get(sid) : cur_swarm_num_different_;
    }

    numSeqs_t SharingSplitFullDacsSwarms::swarm_num_singletons(numSeqs_t sid) const {
        return (finalised_) ? swarms_num_singletons_.get(sid) : cur_swarm_num_singletons_;
    }

    numSeqs_t SharingSplitFullDacsSwarms::swarm_mass(numSeqs_t sid) const {
        return (finalised_) ? swarms_mass_.get(sid) : cur_swarm_mass_;
    }


    /* === SharingSplitFullDacsSwarms::ForwardingDacsSwarm === */

    SharingSplitFullDacsSwarms::ForwardingDacsSwarm::ForwardingDacsSwarm(SharingSplitFullDacsSwarms* swarms, numSeqs_t sid) {

        swarms_ = swarms;
        sid_ = sid;

    }

    SharingSplitFullDacsSwarms::ForwardingDacsSwarm* SharingSplitFullDacsSwarms::ForwardingDacsSwarm::clone() const {
        return new ForwardingDacsSwarm(*this);
    }


    numSeqs_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::seed() const {
        return swarms_->member(sid_, 0);
    }
    SwarmEntry* SharingSplitFullDacsSwarms::ForwardingDacsSwarm::seed_entry() const {
        return new ForwardingSwarmEntry<ForwardingDacsSwarm>(this, 0);
    }
    SwarmEntry* SharingSplitFullDacsSwarms::ForwardingDacsSwarm::seed_entry(SwarmEntry* entry) const {

        dynamic_cast<ForwardingSwarmEntry<ForwardingDacsSwarm>*>(entry)->reassign(this, 0);
        return entry;

    }

    numSeqs_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::get(numSeqs_t i) {
        return swarms_->member(sid_, i);
    }
    SwarmEntry* SharingSplitFullDacsSwarms::ForwardingDacsSwarm::get_entry(numSeqs_t i) {
        return new ForwardingSwarmEntry<ForwardingDacsSwarm>(this, i);
    }
    SwarmEntry* SharingSplitFullDacsSwarms::ForwardingDacsSwarm::get_entry(numSeqs_t i, SwarmEntry* entry) const {

        dynamic_cast<ForwardingSwarmEntry<ForwardingDacsSwarm>*>(entry)->reassign(this, i);
        return entry;

    }

    numSeqs_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::member(numSeqs_t i) const {
        return swarms_->member(sid_, i);
    }

    numSeqs_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::parent(numSeqs_t i) const {
        return swarms_->parent(sid_, i);
    }

    numSeqs_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::gen(numSeqs_t i) const {
        return swarms_->gen(sid_, i);
    }

    dist_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::rad(numSeqs_t i) const {
        return swarms_->rad(sid_, i);
    }

    dist_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::parent_dist(numSeqs_t i) const {
        return swarms_->parent_dist(sid_, i);
    }

    numSeqs_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::size() const {
        return swarms_->swarm_size(sid_);
    }

    numSeqs_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::num_different() const {
        return swarms_->swarm_num_different(sid_);
    }

    numSeqs_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::num_singletons() const {
        return swarms_->swarm_num_singletons(sid_);
    }

    numSeqs_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::mass() const {
        return swarms_->swarm_mass(sid_);
    }

    void SharingSplitFullDacsSwarms::ForwardingDacsSwarm::sort_next_generation(const AmpliconCollection& ac, numSeqs_t p) {
        swarms_->sort_next_generation(ac);
    }

    void SharingSplitFullDacsSwarms::ForwardingDacsSwarm::append(numSeqs_t m, numSeqs_t g, dist_t r, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) {
        swarms_->append(m, p, dist, g, r, ab, new_seq);
    }

    void SharingSplitFullDacsSwarms::ForwardingDacsSwarm::append(numSeqs_t m, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) {

        numSeqs_t cur_pos = swarms_->cur_swarm_start_ + p - swarms_->total_num_members_;
        swarms_->append(m, swarms_->cur_members_[cur_pos], dist, swarms_->cur_gen_ + 1, swarms_->cur_rads_[cur_pos] + dist, ab, new_seq);

    }

    void SharingSplitFullDacsSwarms::ForwardingDacsSwarm::attach(const numSeqs_t p, const numSeqs_t c, Swarm* cs, const dist_t d) {

        if (swarms_->swarms_graftings_[sid_] == nullptr) {
            swarms_->swarms_graftings_[sid_] = new SimpleGraftingInfo();
        }

        swarms_->swarms_graftings_[sid_]->add_link(p, c, cs, d);
        cs->mark_as_attached();

    }

    void SharingSplitFullDacsSwarms::ForwardingDacsSwarm::mark_as_attached() {
        swarms_->swarms_attached_[sid_] = true;
    }

    bool SharingSplitFullDacsSwarms::ForwardingDacsSwarm::is_attached() const {
        return swarms_->swarms_attached_[sid_];
    }

    numSeqs_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::total_size() const {

        numSeqs_t sum = size();
        if (swarms_->swarms_graftings_[sid_] != nullptr) {

            for (numSeqs_t i = 0; i < swarms_->swarms_graftings_[sid_]->num_links(); i++) {
                sum += swarms_->swarms_graftings_[sid_]->get_child_swarm(i)->size();
            }

        }

        return sum;

    }

    numSeqs_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::total_num_different() const {

        numSeqs_t sum = num_different();
        if (swarms_->swarms_graftings_[sid_] != nullptr) {

            for (numSeqs_t i = 0; i < swarms_->swarms_graftings_[sid_]->num_links(); i++) {
                sum += swarms_->swarms_graftings_[sid_]->get_child_swarm(i)->num_different();
            }

        }

        return sum;

    }

    numSeqs_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::total_num_singletons() const {

        numSeqs_t sum = num_singletons();
        if (swarms_->swarms_graftings_[sid_] != nullptr) {

            for (numSeqs_t i = 0; i < swarms_->swarms_graftings_[sid_]->num_links(); i++) {
                sum += swarms_->swarms_graftings_[sid_]->get_child_swarm(i)->num_singletons();
            }

        }

        return sum;

    }

    numSeqs_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::total_mass() const {

        numSeqs_t sum = mass();
        if (swarms_->swarms_graftings_[sid_] != nullptr) {

            for (numSeqs_t i = 0; i < swarms_->swarms_graftings_[sid_]->num_links(); i++) {
                sum += swarms_->swarms_graftings_[sid_]->get_child_swarm(i)->mass();
            }

        }

        return sum;

    }

    dist_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::max_rad() const {
        return swarms_->swarms_max_rad_.get(sid_);
    }

    lenSeqs_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::max_gen() const {
        return (size() > 0) ? swarms_->gen(sid_, size() - 1) : 0;
    }

    GraftingInfo* SharingSplitFullDacsSwarms::ForwardingDacsSwarm::get_grafting_info() const {
        return swarms_->swarms_graftings_[sid_];
    }

    size_t SharingSplitFullDacsSwarms::ForwardingDacsSwarm::size_in_bytes() const {
        return sizeof(ForwardingDacsSwarm); // ForwardingDacsSwarm itself (no ownership over swarms_)
    }

}

