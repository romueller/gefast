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

#include "../include/Swarms.hpp"
#include "../include/space/Basics.hpp"

namespace GeFaST {

    /* === SimpleSwarm === */

    SimpleSwarm::SimpleSwarm(numSeqs_t m, numSeqs_t ab) {

        members_.push_back(m);
        parents_.push_back(m);
        parent_dists_.push_back(0);
        gens_.push_back(0);
        rads_.push_back(0);

        mass_ = ab;
        num_different_ = 1;
        num_singletons_ = (ab == 1);

        max_rad_ = 0;

        graftings_ = nullptr;
        attached_ = false;

    }

    SimpleSwarm::~SimpleSwarm() {
        delete graftings_;
    }

    SimpleSwarm::SimpleSwarm(const SimpleSwarm& other) : // copy constructor
            members_(other.members_), parents_(other.parents_), parent_dists_(other.parent_dists_),
            gens_(other.gens_), rads_(other.rads_), mass_(other.mass_), num_different_(other.num_different_),
            num_singletons_(other.num_singletons_), max_rad_(other.max_rad_),
            graftings_(other.graftings_ ? other.graftings_->clone() : nullptr), attached_(other.attached_) {

        // nothing else to do

    }

    SimpleSwarm::SimpleSwarm(SimpleSwarm&& other) noexcept : // move constructor
            members_(std::move(other.members_)), parents_(std::move(other.parents_)),
            parent_dists_(std::move(other.parent_dists_)), gens_(std::move(other.gens_)), rads_(std::move(other.rads_)),
            mass_(other.mass_), num_different_(other.num_different_), num_singletons_(other.num_singletons_),
            max_rad_(other.max_rad_), graftings_(other.graftings_), attached_(other.attached_) {

        other.graftings_ = nullptr;

    }

    SimpleSwarm& SimpleSwarm::operator=(const SimpleSwarm& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete graftings_;

        // copy new resources
        members_ = other.members_;
        parents_ = other.parents_;
        parent_dists_ = other.parent_dists_;
        gens_ = other.gens_;
        rads_ = other.rads_;

        mass_ = other.mass_;
        num_different_ = other.num_different_;
        num_singletons_ = other.num_singletons_;

        max_rad_ = other.max_rad_;

        graftings_ = other.graftings_ ? other.graftings_->clone() : nullptr;
        attached_ = other.attached_;

        return *this;

    }

    SimpleSwarm& SimpleSwarm::operator=(SimpleSwarm&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete graftings_;

        // copy new resources
        members_ = other.members_;
        parents_ = other.parents_;
        parent_dists_ = other.parent_dists_;
        gens_ = other.gens_;
        rads_ = other.rads_;

        mass_ = other.mass_;
        num_different_ = other.num_different_;
        num_singletons_ = other.num_singletons_;

        max_rad_ = other.max_rad_;

        graftings_ = other.graftings_; other.graftings_ = nullptr;
        attached_ = other.attached_;

        return *this;

    }

    SimpleSwarm* SimpleSwarm::clone() const {
        return new SimpleSwarm(*this);
    }

    numSeqs_t SimpleSwarm::seed() const {
        return members_[0];
    }

    SwarmEntry* SimpleSwarm::seed_entry() const {
        return new ForwardingSwarmEntry<SimpleSwarm>(this, 0);
    }

    SwarmEntry* SimpleSwarm::seed_entry(SwarmEntry* entry) const {

        dynamic_cast<ForwardingSwarmEntry<SimpleSwarm>*>(entry)->reassign(this, 0);
        return entry;

    }

    numSeqs_t SimpleSwarm::get(numSeqs_t i) {
        return members_[i];
    }

    SwarmEntry* SimpleSwarm::get_entry(numSeqs_t i) {
        return new ForwardingSwarmEntry<SimpleSwarm>(this, i);
    }

    SwarmEntry* SimpleSwarm::get_entry(numSeqs_t i, SwarmEntry* entry) const {

        dynamic_cast<ForwardingSwarmEntry<SimpleSwarm>*>(entry)->reassign(this, i);
        return entry;

    }

    numSeqs_t SimpleSwarm::member(numSeqs_t i) const {
        return members_[i];
    }

    numSeqs_t SimpleSwarm::parent(numSeqs_t i) const {
        return parents_[i];
    }

    numSeqs_t SimpleSwarm::gen(numSeqs_t i) const {
        return gens_[i];
    }

    dist_t SimpleSwarm::rad(numSeqs_t i) const {
        return rads_[i];
    }

    dist_t SimpleSwarm::parent_dist(numSeqs_t i) const {
        return parent_dists_[i];
    }

    numSeqs_t SimpleSwarm::size() const {
        return members_.size();
    }

    numSeqs_t SimpleSwarm::num_different() const {
        return num_different_;
    }

    numSeqs_t SimpleSwarm::num_singletons() const {
        return num_singletons_;
    }

    numSeqs_t SimpleSwarm::mass() const {
        return mass_;
    }

    void SimpleSwarm::sort_next_generation(const AmpliconCollection& ac, numSeqs_t p) {

        std::vector<numSeqs_t> permutation(size() - p);
        std::iota(permutation.begin(), permutation.end(), p);

        std::sort(permutation.begin(), permutation.end(),
                  [&](const numSeqs_t a, const numSeqs_t b) {
                      return (ac.ab(members_[a]) > ac.ab(members_[b])) ||
                             ((ac.ab(members_[a]) == ac.ab(members_[b])) &&
                              (*ac.id(members_[a]) < *ac.id(members_[b])));
                  }
        );

        std::vector<bool> done(size());
        for (numSeqs_t i = p; i < size(); ++i) {

            if (!done[i]) {

                done[i] = true;
                numSeqs_t prev_j = i;
                numSeqs_t j = permutation[i - p];

                while (i != j) {

                    std::swap(members_[prev_j], members_[j]);
                    std::swap(parents_[prev_j], parents_[j]);
                    std::swap(parent_dists_[prev_j], parent_dists_[j]);
                    // no swaps necessary for gens_
                    std::swap(rads_[prev_j], rads_[j]);
                    done[j] = true;
                    prev_j = j;
                    j = permutation[j - p];

                }
            }

        }

    }

    void SimpleSwarm::append(numSeqs_t m, numSeqs_t g, dist_t r, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) {

        members_.push_back(m);
        parents_.push_back(p);
        parent_dists_.push_back(dist);
        gens_.push_back(g);
        rads_.push_back(r);

        mass_ += ab;
        num_different_ += new_seq;
        num_singletons_ += (ab == 1);

        max_rad_ = std::max(max_rad_, rads_.back());

    }

    void SimpleSwarm::append(numSeqs_t m, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) {

        members_.push_back(m);
        parents_.push_back(members_[p]);
        parent_dists_.push_back(dist);
        gens_.push_back(gens_[p] + 1);
        rads_.push_back(rads_[p] + dist);

        mass_ += ab;
        num_different_ += new_seq;
        num_singletons_ += (ab == 1);

        max_rad_ = std::max(max_rad_, rads_.back());

    }

    void SimpleSwarm::attach(const numSeqs_t p, const numSeqs_t c, Swarm* cs, const dist_t d) {

        if (graftings_ == nullptr) {
            graftings_ = new SimpleGraftingInfo();
        }

        graftings_->add_link(p, c, cs, d);
        cs->mark_as_attached();

    }

    void SimpleSwarm::mark_as_attached() {
        attached_ = true;
    }

    bool SimpleSwarm::is_attached() const {
        return attached_;
    }

    numSeqs_t SimpleSwarm::total_size() const {

        numSeqs_t sum = size();
        if (graftings_ != nullptr) {

            for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                sum += graftings_->get_child_swarm(i)->size();
            }

        }

        return sum;

    }

    numSeqs_t SimpleSwarm::total_num_different() const {

        numSeqs_t sum = num_different();
        if (graftings_ != nullptr) {

            for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                sum += graftings_->get_child_swarm(i)->num_different();
            }

        }

        return sum;

    }

    numSeqs_t SimpleSwarm::total_num_singletons() const {

        numSeqs_t sum = num_singletons();
        if (graftings_ != nullptr) {

            for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                sum += graftings_->get_child_swarm(i)->num_singletons();
            }

        }

        return sum;

    }

    numSeqs_t SimpleSwarm::total_mass() const {

        numSeqs_t sum = mass();
        if (graftings_ != nullptr) {

            for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                sum += graftings_->get_child_swarm(i)->mass();
            }

        }

        return sum;

    }

    dist_t SimpleSwarm::max_rad() const {
        return max_rad_;
    }

    lenSeqs_t SimpleSwarm::max_gen() const {
        return gens_.back();
    }

    GraftingInfo* SimpleSwarm::get_grafting_info() const {
        return graftings_;
    }

    size_t SimpleSwarm::size_in_bytes() const {

        std::cerr << "ERROR: SimpleSwarm is currently not part of the space-efficiency analysis." << std::endl;
        return 0;

    }


    /* === SimpleCombinedSwarm === */

    SimpleCombinedSwarm::SimpleCombinedSwarm(numSeqs_t m, numSeqs_t ab) {

        members_.emplace_back(m, m, 0, 0, 0);

        mass_ = ab;
        num_different_ = 1;
        num_singletons_ = (ab == 1);

        max_rad_ = 0;

        graftings_ = nullptr;
        attached_ = false;

    }

    SimpleCombinedSwarm::~SimpleCombinedSwarm() {
        delete graftings_;
    }

    SimpleCombinedSwarm::SimpleCombinedSwarm(const SimpleCombinedSwarm& other) : // copy constructor
            members_(other.members_), mass_(other.mass_), num_different_(other.num_different_),
            num_singletons_(other.num_singletons_), max_rad_(other.max_rad_),
            graftings_(other.graftings_ ? other.graftings_->clone() : nullptr), attached_(other.attached_) {

        // nothing else to do

    }

    SimpleCombinedSwarm::SimpleCombinedSwarm(SimpleCombinedSwarm&& other) noexcept : // move constructor
            members_(std::move(other.members_)), mass_(other.mass_), num_different_(other.num_different_),
            num_singletons_(other.num_singletons_), max_rad_(other.max_rad_), graftings_(other.graftings_), attached_(other.attached_) {

        other.graftings_ = nullptr;

    }

    SimpleCombinedSwarm& SimpleCombinedSwarm::operator=(const SimpleCombinedSwarm& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete graftings_;

        // copy new resources
        members_ = other.members_;

        mass_ = other.mass_;
        num_different_ = other.num_different_;
        num_singletons_ = other.num_singletons_;

        max_rad_ = other.max_rad_;

        graftings_ = other.graftings_ ? other.graftings_->clone() : nullptr;
        attached_ = other.attached_;

        return *this;

    }

    SimpleCombinedSwarm& SimpleCombinedSwarm::operator=(SimpleCombinedSwarm&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete graftings_;

        // copy new resources
        members_ = other.members_;

        mass_ = other.mass_;
        num_different_ = other.num_different_;
        num_singletons_ = other.num_singletons_;

        max_rad_ = other.max_rad_;

        graftings_ = other.graftings_; other.graftings_ = nullptr;
        attached_ = other.attached_;

        return *this;

    }

    SimpleCombinedSwarm* SimpleCombinedSwarm::clone() const {
        return new SimpleCombinedSwarm(*this);
    }

    numSeqs_t SimpleCombinedSwarm::seed() const {
        return members_[0].member;
    }

    SwarmEntry* SimpleCombinedSwarm::seed_entry() const {
        return new ForwardingSwarmEntry<SimpleCombinedSwarm>(this, 0);
    }

    SwarmEntry* SimpleCombinedSwarm::seed_entry(SwarmEntry* entry) const {

        dynamic_cast<ForwardingSwarmEntry<SimpleCombinedSwarm>*>(entry)->reassign(this, 0);
        return entry;

    }

    numSeqs_t SimpleCombinedSwarm::get(numSeqs_t i) {
        return members_[i].member;
    }

    SwarmEntry* SimpleCombinedSwarm::get_entry(numSeqs_t i) {
        return new ForwardingSwarmEntry<SimpleCombinedSwarm>(this, i);
    }

    SwarmEntry* SimpleCombinedSwarm::get_entry(numSeqs_t i, SwarmEntry* entry) const {

        dynamic_cast<ForwardingSwarmEntry<SimpleCombinedSwarm>*>(entry)->reassign(this, i);
        return entry;

    }

    numSeqs_t SimpleCombinedSwarm::member(numSeqs_t i) const {
        return members_[i].member;
    }

    numSeqs_t SimpleCombinedSwarm::parent(numSeqs_t i) const {
        return members_[i].parent;
    }

    numSeqs_t SimpleCombinedSwarm::gen(numSeqs_t i) const {
        return members_[i].gen;
    }

    dist_t SimpleCombinedSwarm::rad(numSeqs_t i) const {
        return members_[i].rad;
    }

    dist_t SimpleCombinedSwarm::parent_dist(numSeqs_t i) const {
        return members_[i].parent_dist;
    }

    numSeqs_t SimpleCombinedSwarm::size() const {
        return members_.size();
    }

    numSeqs_t SimpleCombinedSwarm::num_different() const {
        return num_different_;
    }

    numSeqs_t SimpleCombinedSwarm::num_singletons() const {
        return num_singletons_;
    }

    numSeqs_t SimpleCombinedSwarm::mass() const {
        return mass_;
    }

    void SimpleCombinedSwarm::sort_next_generation(const AmpliconCollection& ac, numSeqs_t p) {

        std::vector<numSeqs_t> permutation(size() - p);
        std::iota(permutation.begin(), permutation.end(), p);

        std::sort(permutation.begin(), permutation.end(),
                  [&](const numSeqs_t a, const numSeqs_t b) {
                      return (ac.ab(members_[a].member) > ac.ab(members_[b].member)) ||
                             ((ac.ab(members_[a].member) == ac.ab(members_[b].member)) &&
                              (*ac.id(members_[a].member) < *ac.id(members_[b].member)));
                  }
        );

        std::vector<bool> done(size());
        for (numSeqs_t i = p; i < size(); ++i) {

            if (!done[i]) {

                done[i] = true;
                numSeqs_t prev_j = i;
                numSeqs_t j = permutation[i - p];

                while (i != j) {

                    std::swap(members_[prev_j], members_[j]);

                    done[j] = true;
                    prev_j = j;
                    j = permutation[j - p];

                }
            }

        }

    }

    void SimpleCombinedSwarm::append(numSeqs_t m, numSeqs_t g, dist_t r, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) {

        members_.emplace_back(m, p, dist, g, r);

        mass_ += ab;
        num_different_ += new_seq;
        num_singletons_ += (ab == 1);

        max_rad_ = std::max(max_rad_, r);

    }

    void SimpleCombinedSwarm::append(numSeqs_t m, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) {

        members_.emplace_back(m, members_[p].member, dist, members_[p].gen + 1, members_[p].rad + dist);

        mass_ += ab;
        num_different_ += new_seq;
        num_singletons_ += (ab == 1);

        max_rad_ = std::max(max_rad_, members_.back().rad);

    }

    void SimpleCombinedSwarm::attach(const numSeqs_t p, const numSeqs_t c, Swarm* cs, const dist_t d) {

        if (graftings_ == nullptr) {
            graftings_ = new SimpleGraftingInfo();
        }

        graftings_->add_link(p, c, cs, d);
        cs->mark_as_attached();

    }

    void SimpleCombinedSwarm::mark_as_attached() {
        attached_ = true;
    }

    bool SimpleCombinedSwarm::is_attached() const {
        return attached_;
    }

    numSeqs_t SimpleCombinedSwarm::total_size() const {

        numSeqs_t sum = size();
        if (graftings_ != nullptr) {

            for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                sum += graftings_->get_child_swarm(i)->size();
            }

        }

        return sum;

    }

    numSeqs_t SimpleCombinedSwarm::total_num_different() const {

        numSeqs_t sum = num_different();
        if (graftings_ != nullptr) {

            for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                sum += graftings_->get_child_swarm(i)->num_different();
            }

        }

        return sum;

    }

    numSeqs_t SimpleCombinedSwarm::total_num_singletons() const {

        numSeqs_t sum = num_singletons();
        if (graftings_ != nullptr) {

            for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                sum += graftings_->get_child_swarm(i)->num_singletons();
            }

        }

        return sum;

    }

    numSeqs_t SimpleCombinedSwarm::total_mass() const {

        numSeqs_t sum = mass();
        if (graftings_ != nullptr) {

            for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                sum += graftings_->get_child_swarm(i)->mass();
            }

        }

        return sum;

    }

    dist_t SimpleCombinedSwarm::max_rad() const {
        return max_rad_;
    }

    lenSeqs_t SimpleCombinedSwarm::max_gen() const {
        return members_.back().gen;
    }

    GraftingInfo* SimpleCombinedSwarm::get_grafting_info() const {
        return graftings_;
    }

    size_t SimpleCombinedSwarm::size_in_bytes() const {

        std::cerr << "ERROR: SimpleCombinedSwarm is currently not part of the space-efficiency analysis." << std::endl;
        return 0;

    }

}