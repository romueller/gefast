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

#ifndef GEFAST_SWARMS_HPP
#define GEFAST_SWARMS_HPP

#include <algorithm>
#include <cstring>

#include "Base.hpp"

namespace GeFaST {

    /*
     * Simple representation of a swarm (cluster).
     *
     * The different information on each member of the swarm (parent, generation etc.)
     * are stored in respective STL vectors.
     * Stores some overall characteristics (e.g. mass) explicitly.
     * Swarms grafted to a particular swarm instance are only added implicitly
     * through information on the grafting links.
     */
    class SimpleSwarm : public Swarm {
    public:
        /*
         * Nested class for conveniently handling individual swarm members and accessing its information.
         */
        class SimpleSwarmEntry : public SwarmEntry {

            friend class SimpleSwarm;

        public:
            SimpleSwarmEntry(const SimpleSwarm* swarm, numSeqs_t i) : swarm_(swarm), internal_id_(i) {
                // nothing else to do
            }

            virtual ~SimpleSwarmEntry() = default;

            SimpleSwarmEntry(const SimpleSwarmEntry& other) = default; // copy constructor

            SimpleSwarmEntry(SimpleSwarmEntry&& other) = default; // move constructor

            SimpleSwarmEntry& operator=(const SimpleSwarmEntry& other) = default; // copy assignment operator

            SimpleSwarmEntry& operator=(SimpleSwarmEntry&& other) = default; // move assignment operator

            SimpleSwarmEntry* clone() const override { // deep-copy clone method
                return new SimpleSwarmEntry(*this);
            }

            inline numSeqs_t member() const override {
                return swarm_->member(internal_id_);
            }

            inline numSeqs_t parent() const override {
                return swarm_->parent(internal_id_);
            }

            inline numSeqs_t gen() const override {
                return swarm_->gen(internal_id_);
            }

            inline dist_t rad() const override {
                return swarm_->rad(internal_id_);
            }

            inline dist_t parent_dist() const override {
                return swarm_->parent_dist(internal_id_);
            }

        protected:
            SimpleSwarmEntry() = default;

            const SimpleSwarm* swarm_; // no ownership
            numSeqs_t internal_id_;

        };

        /*
         * Create new swarm and initialise it with member m of abundance ab as the seed.
         */
        SimpleSwarm(numSeqs_t m, numSeqs_t ab);

        virtual ~SimpleSwarm();

        SimpleSwarm(const SimpleSwarm& other); // copy constructor

        SimpleSwarm(SimpleSwarm&& other) noexcept; // move constructor

        SimpleSwarm& operator=(const SimpleSwarm& other); // copy assignment operator

        SimpleSwarm& operator=(SimpleSwarm&& other) noexcept; // move assignment operator

        SimpleSwarm* clone() const override; // deep-copy clone method


        numSeqs_t seed() const override;
        SwarmEntry* seed_entry() const override;
        SwarmEntry* seed_entry(SwarmEntry* entry) const override;

        numSeqs_t get(numSeqs_t i) override;
        SwarmEntry* get_entry(numSeqs_t i) override;
        SwarmEntry* get_entry(numSeqs_t i, SwarmEntry* entry) const override;

        numSeqs_t member(numSeqs_t i) const override;

        numSeqs_t parent(numSeqs_t i) const override;

        numSeqs_t gen(numSeqs_t i) const override;

        dist_t rad(numSeqs_t i) const override;

        dist_t parent_dist(numSeqs_t i) const override;

        numSeqs_t size() const override;

        numSeqs_t num_different() const override;

        numSeqs_t num_singletons() const override;

        numSeqs_t mass() const override;

        void sort_next_generation(const AmpliconCollection& ac, numSeqs_t p) override;

        void append(numSeqs_t m, numSeqs_t g, dist_t r, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) override;

        void append(numSeqs_t m, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) override;

        void attach(const numSeqs_t p, const numSeqs_t c, Swarm* cs, const dist_t d) override;

        void mark_as_attached() override;

        bool is_attached() const override;

        numSeqs_t total_size() const override;

        numSeqs_t total_num_different() const override;

        numSeqs_t total_num_singletons() const override;

        numSeqs_t total_mass() const override;

        dist_t max_rad() const override;

        lenSeqs_t max_gen() const override;

        GraftingInfo* get_grafting_info() const override;

    protected:
        SimpleSwarm() = default;

        // the information on each swarm member are stored separately but in the same order in multiple containers
        std::vector<numSeqs_t> members_; // pool-internal integer id of swarm member
        std::vector<numSeqs_t> parents_; // pool-internal integer id of its parent
        std::vector<dist_t> parent_dists_; // distance to its parent
        std::vector<numSeqs_t> gens_; // its generation number
        std::vector<dist_t> rads_; // its radius (accumulated distance to root)

        // overall characteristics of the swarm (not considering the grafted amplicons)
        numSeqs_t mass_; // sum of the abundances of amplicons
        numSeqs_t num_different_; // number of different sequences
        numSeqs_t num_singletons_; // number of singletons
        dist_t max_rad_; // maximum radius

        GraftingInfo* graftings_; // grafting links with this swarm as the parent (nullptr when there are none)
        bool attached_; // Boolean flag indicating whether this swarm is attached (as a child) to another swarm

    };


    /*
     * Simple representation of a collection of swarms.
     *
     * The swarms are simply stored as separate instances in an STL vector.
     */
    template<typename S>
    class SimpleSwarms : public Swarms {

    public:
        SimpleSwarms() = default;

        Swarm& initialise_cluster(numSeqs_t s, numSeqs_t ab) override {

            swarms_.push_back(new S(s, ab));
            return *(swarms_.back());

        }

        ~SimpleSwarms() {

            for (auto i = 0; i < swarms_.size(); i++) {
                delete swarms_[i]; swarms_[i] = nullptr;
            }

        }

        SimpleSwarms(const SimpleSwarms& other) { // copy constructor

            for (auto s : other.swarms_) {
                swarms_.push_back(s->clone());
            }

        }

        SimpleSwarms(SimpleSwarms&& other) : swarms_(std::move(other.swarms_)) { // move constructor
            // nothing else to do
        }

        SimpleSwarms& operator=(const SimpleSwarms& other) {// copy assignment operator

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

        }

        SimpleSwarms& operator=(SimpleSwarms&& other) { // move assignment operator

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

        }

        SimpleSwarms* clone() const override { // deep-copy clone method
            return new SimpleSwarms(*this);
        }


        numSeqs_t size() const override {
            return swarms_.size();
        }

        Swarm& get(const numSeqs_t i) override {
            return *(swarms_[i]);
        }

        const Swarm& get(const numSeqs_t i) const override {
            return *(swarms_[i]);
        }

        void clear() override {

            for (auto i = 0; i < swarms_.size(); i++) {
                delete swarms_[i]; swarms_[i] = nullptr;
            }
            swarms_.clear();

        }

    private:
        std::vector<Swarm*> swarms_; // swarms contained in the collection

    };

}

#endif //GEFAST_SWARMS_HPP
