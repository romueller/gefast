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

#ifndef GEFAST_SWARMS_HPP
#define GEFAST_SWARMS_HPP

#include <algorithm>
#include <cstring>
#include <numeric>

#include "Base.hpp"

namespace GeFaST {

    /*
     * Class for conveniently handling individual swarm members and accessing its information without storing a copy.
     */
    template<typename S>
    class ForwardingSwarmEntry : public SwarmEntry {

    public:
        ForwardingSwarmEntry(const S* swarm, numSeqs_t i) : swarm_(swarm), internal_id_(i) {
            // nothing else to do
        }

        virtual ~ForwardingSwarmEntry() = default;

        ForwardingSwarmEntry(const ForwardingSwarmEntry& other) = default; // copy constructor

        ForwardingSwarmEntry(ForwardingSwarmEntry&& other) = default; // move constructor

        ForwardingSwarmEntry& operator=(const ForwardingSwarmEntry& other) = default; // copy assignment operator

        ForwardingSwarmEntry& operator=(ForwardingSwarmEntry&& other) = default; // move assignment operator

        ForwardingSwarmEntry* clone() const override { // deep-copy clone method
            return new ForwardingSwarmEntry(*this);
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

        inline void reassign(const S* s, numSeqs_t i) {

            swarm_ = s;
            internal_id_ = i;

        }

    protected:
        ForwardingSwarmEntry() = default;

        const S* swarm_; // no ownership
        numSeqs_t internal_id_;

    };



    /* ==== Swarms implementations ==== */

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


    /*
     * Representation of a collection of swarms.
     *
     * The swarms and their overall characteristics (e.g. mass) are stored as separate instances
     * but the members of the different swarms are held collectively in the Swarms instance.
     * The different Swarm instances just point to their members.
     */
    class SharingSwarms : public Swarms {

    public:
        virtual ~SharingSwarms() = default;

        virtual numSeqs_t num_members() const = 0;

        virtual void append(numSeqs_t m, numSeqs_t p, dist_t d, numSeqs_t g, dist_t r) = 0;

        virtual numSeqs_t get_member(numSeqs_t i) const = 0;

        virtual numSeqs_t get_parent(numSeqs_t i) const = 0;

        virtual dist_t get_parent_dist(numSeqs_t i) const = 0;

        virtual numSeqs_t get_gen(numSeqs_t i) const = 0;

        virtual dist_t get_rad(numSeqs_t i) const = 0;

        virtual void swap(numSeqs_t i, numSeqs_t j) = 0;

    };

    /*
     * The information of each swarm member are split and stored separately (but in the same order)
     * in different vectors.
     */
    template<typename S>
    class SharingSplitSwarms : public SharingSwarms {

    public:
        SharingSplitSwarms() = default;

        SharingSplitSwarms(const numSeqs_t pool_size) {

            all_members_.reserve(pool_size);
            all_parents_.reserve(pool_size);
            all_parent_dists_.reserve(pool_size);
            all_gens_.reserve(pool_size);
            all_rads_.reserve(pool_size);

        }

        Swarm& initialise_cluster(numSeqs_t s, numSeqs_t ab) override {

            swarms_.push_back(new S(*this, s, ab));
            return *(swarms_.back());

        }

        ~SharingSplitSwarms() {

            for (auto i = 0; i < swarms_.size(); i++) {
                delete swarms_[i]; swarms_[i] = nullptr;
            }

        }

        SharingSplitSwarms(const SharingSplitSwarms& other) : all_members_(other.all_members_), all_parents_(other.all_parents_),
                all_parent_dists_(other.all_parent_dists_), all_gens_(other.all_gens_), all_rads_(other.all_rads_) { // copy constructor

            swarms_.reserve(other.size());
            for (auto s : other.swarms_) {
                swarms_.push_back(s->clone());
            }

        }

        SharingSplitSwarms(SharingSplitSwarms&& other) : swarms_(std::move(other.swarms_)), all_members_(std::move(other.all_members_)),
                all_parents_(std::move(other.all_parents_)), all_parent_dists_(std::move(other.all_parent_dists_)),
                all_gens_(std::move(other.all_gens_)), all_rads_(std::move(other.all_rads_)) { // move constructor

            // nothing else to do

        }

        SharingSplitSwarms& operator=(const SharingSplitSwarms& other) {// copy assignment operator

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

        }

        SharingSplitSwarms& operator=(SharingSplitSwarms&& other) { // move assignment operator

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

        }

        SharingSplitSwarms* clone() const override { // deep-copy clone method
            return new SharingSplitSwarms(*this);
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

            all_members_.clear();
            all_parents_.clear();
            all_parent_dists_.clear();
            all_gens_.clear();
            all_rads_.clear();

        }

        numSeqs_t num_members() const override {
            return all_members_.size();
        }

        void append(numSeqs_t m, numSeqs_t p, dist_t d, numSeqs_t g, dist_t r) override {

            all_members_.emplace_back(m);
            all_parents_.emplace_back(p);
            all_parent_dists_.emplace_back(d);
            all_gens_.emplace_back(g);
            all_rads_.emplace_back(r);

        }

        inline numSeqs_t get_member(numSeqs_t i) const override {
            return all_members_[i];
        }

        inline numSeqs_t get_parent(numSeqs_t i) const override {
            return all_parents_[i];
        }

        inline dist_t get_parent_dist(numSeqs_t i) const override {
            return all_parent_dists_[i];
        }

        inline numSeqs_t get_gen(numSeqs_t i) const override {
            return all_gens_[i];
        }

        inline dist_t get_rad(numSeqs_t i) const override {
            return all_rads_[i];
        }

        inline void swap(numSeqs_t i, numSeqs_t j) override {

            std::swap(all_members_[i], all_members_[j]);
            std::swap(all_parents_[i], all_parents_[j]);
            std::swap(all_parent_dists_[i], all_parent_dists_[j]);
            // no swaps necessary for gens_
            std::swap(all_rads_[i], all_rads_[j]);

        }

    private:
        std::vector<Swarm*> swarms_; // swarms contained in the collection

        std::vector<numSeqs_t> all_members_;
        std::vector<numSeqs_t> all_parents_;
        std::vector<dist_t> all_parent_dists_;
        std::vector<numSeqs_t> all_gens_;
        std::vector<dist_t> all_rads_;

    };

    /*
     * The information of each swarm member are stored together as a Member struct in a single vector.
     */
    template<typename S>
    class SharingCombinedSwarms : public SharingSwarms {

        struct Member {

            Member() = default;

            Member(numSeqs_t m, numSeqs_t p, dist_t d, numSeqs_t g, dist_t r)
                    : member(m), parent(p), parent_dist(d), gen(g), rad(r) {

                // nothing else to do

            }

            numSeqs_t member;
            numSeqs_t parent;
            dist_t parent_dist;
            numSeqs_t gen;
            dist_t rad;

        };

    public:
        SharingCombinedSwarms() = default;

        SharingCombinedSwarms(const numSeqs_t pool_size) {
            all_members_.reserve(pool_size);
        }

        Swarm& initialise_cluster(numSeqs_t s, numSeqs_t ab) override {

            swarms_.push_back(new S(*this, s, ab));
            return *(swarms_.back());

        }

        ~SharingCombinedSwarms() {

            for (auto i = 0; i < swarms_.size(); i++) {
                delete swarms_[i]; swarms_[i] = nullptr;
            }

        }

        SharingCombinedSwarms(const SharingCombinedSwarms& other) : all_members_(other.all_members_) { // copy constructor

            swarms_.reserve(other.size());
            for (auto s : other.swarms_) {
                swarms_.push_back(s->clone());
            }

        }

        SharingCombinedSwarms(SharingCombinedSwarms&& other) : swarms_(std::move(other.swarms_)), all_members_(std::move(other.all_members_)) { // move constructor

            // nothing else to do

        }

        SharingCombinedSwarms& operator=(const SharingCombinedSwarms& other) {// copy assignment operator

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

        }

        SharingCombinedSwarms& operator=(SharingCombinedSwarms&& other) { // move assignment operator

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

        }

        SharingCombinedSwarms* clone() const override { // deep-copy clone method
            return new SharingCombinedSwarms(*this);
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

            all_members_.clear();

        }

        inline numSeqs_t num_members() const override {
            return all_members_.size();
        }

        inline void append(numSeqs_t m, numSeqs_t p, dist_t d, numSeqs_t g, dist_t r) {
            all_members_.emplace_back(m, p, d, g, r);
        }

        inline numSeqs_t get_member(numSeqs_t i) const {
            return all_members_[i].member;
        }

        inline numSeqs_t get_parent(numSeqs_t i) const {
            return all_members_[i].parent;
        }

        inline dist_t get_parent_dist(numSeqs_t i) const {
            return all_members_[i].parent_dist;
        }

        inline numSeqs_t get_gen(numSeqs_t i) const {
            return all_members_[i].gen;
        }

        inline dist_t get_rad(numSeqs_t i) const {
            return all_members_[i].rad;
        }

        inline void swap(numSeqs_t i, numSeqs_t j) {
            std::swap(all_members_[i], all_members_[j]);
        }

    private:
        std::vector<Swarm*> swarms_; // swarms contained in the collection

        std::vector<Member> all_members_;

    };



    /* ==== Swarm implementations ==== */

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
         * Create new swarm and initialise it with member m of abundance ab as the seed.
         */
        SimpleSwarm(numSeqs_t m, numSeqs_t ab);

        virtual ~SimpleSwarm();

        SimpleSwarm(const SimpleSwarm& other); // copy constructor

        SimpleSwarm(SimpleSwarm&& other) noexcept; // move constructor

        SimpleSwarm& operator=(const SimpleSwarm& other); // copy assignment operator

        SimpleSwarm& operator=(SimpleSwarm&& other) noexcept; // move assignment operator

        SimpleSwarm* clone() const override; // deep-copy clone method


        inline numSeqs_t seed() const override;
        inline SwarmEntry* seed_entry() const override;
        inline SwarmEntry* seed_entry(SwarmEntry* entry) const override;

        inline numSeqs_t get(numSeqs_t i) override;
        inline SwarmEntry* get_entry(numSeqs_t i) override;
        inline SwarmEntry* get_entry(numSeqs_t i, SwarmEntry* entry) const override;

        inline numSeqs_t member(numSeqs_t i) const override;

        inline numSeqs_t parent(numSeqs_t i) const override;

        inline numSeqs_t gen(numSeqs_t i) const override;

        inline dist_t rad(numSeqs_t i) const override;

        inline dist_t parent_dist(numSeqs_t i) const override;

        inline numSeqs_t size() const override;

        inline numSeqs_t num_different() const override;

        inline numSeqs_t num_singletons() const override;

        inline numSeqs_t mass() const override;

        void sort_next_generation(const AmpliconCollection& ac, numSeqs_t p) override;

        inline void append(numSeqs_t m, numSeqs_t g, dist_t r, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) override;

        inline void append(numSeqs_t m, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) override;

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
     * Simple representation of a swarm (cluster).
     *
     * The different information on each member of the swarm (parent, generation etc.)
     * are stored together as a Member struct in a single vector.
     * Stores some overall characteristics (e.g. mass) explicitly.
     * Swarms grafted to a particular swarm instance are only added implicitly
     * through information on the grafting links.
     */
    class SimpleCombinedSwarm : public Swarm {

        struct Member {

            Member() = default;

            Member(numSeqs_t m, numSeqs_t p, dist_t d, numSeqs_t g, dist_t r)
                    : member(m), parent(p), parent_dist(d), gen(g), rad(r) {

                // nothing else to do

            }

            numSeqs_t member; // pool-internal integer id of swarm member
            numSeqs_t parent; // pool-internal integer id of its parent
            dist_t parent_dist; // distance to its parent
            numSeqs_t gen; // its generation number
            dist_t rad; // its radius (accumulated distance to root)

        };

    public:
        /*
         * Create new swarm and initialise it with member m of abundance ab as the seed.
         */
        SimpleCombinedSwarm(numSeqs_t m, numSeqs_t ab);

        virtual ~SimpleCombinedSwarm();

        SimpleCombinedSwarm(const SimpleCombinedSwarm& other); // copy constructor

        SimpleCombinedSwarm(SimpleCombinedSwarm&& other) noexcept; // move constructor

        SimpleCombinedSwarm& operator=(const SimpleCombinedSwarm& other); // copy assignment operator

        SimpleCombinedSwarm& operator=(SimpleCombinedSwarm&& other) noexcept; // move assignment operator

        SimpleCombinedSwarm* clone() const override; // deep-copy clone method


        inline numSeqs_t seed() const override;
        inline SwarmEntry* seed_entry() const override;
        inline SwarmEntry* seed_entry(SwarmEntry* entry) const override;

        inline numSeqs_t get(numSeqs_t i) override;
        inline SwarmEntry* get_entry(numSeqs_t i) override;
        inline SwarmEntry* get_entry(numSeqs_t i, SwarmEntry* entry) const override;

        inline numSeqs_t member(numSeqs_t i) const override;

        inline numSeqs_t parent(numSeqs_t i) const override;

        inline numSeqs_t gen(numSeqs_t i) const override;

        inline dist_t rad(numSeqs_t i) const override;

        inline dist_t parent_dist(numSeqs_t i) const override;

        inline numSeqs_t size() const override;

        inline numSeqs_t num_different() const override;

        inline numSeqs_t num_singletons() const override;

        inline numSeqs_t mass() const override;

        void sort_next_generation(const AmpliconCollection& ac, numSeqs_t p) override;

        inline void append(numSeqs_t m, numSeqs_t g, dist_t r, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) override;

        inline void append(numSeqs_t m, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) override;

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
        SimpleCombinedSwarm() = default;

        std::vector<Member> members_; // the information on each swarm member are stored collectively

        // overall characteristics of the swarm (not considering the grafted amplicons)
        numSeqs_t mass_; // sum of the abundances of amplicons
        numSeqs_t num_different_; // number of different sequences
        numSeqs_t num_singletons_; // number of singletons
        dist_t max_rad_; // maximum radius

        GraftingInfo* graftings_; // grafting links with this swarm as the parent (nullptr when there are none)
        bool attached_; // Boolean flag indicating whether this swarm is attached (as a child) to another swarm

    };

    /*
     * Representation of a swarm storing only the overall characteristics (e.g. mass).
     * The information on each member of the swarm (parent, generation etc.) are stored in SharingSwarms instance.
     * Requests of such information are forwarded to the SharingSwarms instance.
     * Swarms grafted to a particular swarm instance are only added implicitly
     * through information on the grafting links.
     */
    class ForwardingSwarm : public Swarm {
    public:

        /*
         * Create new swarm and initialise it with member m of abundance ab as the seed.
         */
        ForwardingSwarm(SharingSwarms& as, numSeqs_t m, numSeqs_t ab) {

            swarms_ = &as;
            start_ = as.num_members();

            as.append(m, m, 0, 0, 0);

            size_ = 1;
            mass_ = ab;
            num_different_ = 1;
            num_singletons_ = (ab == 1);

            max_rad_ = 0;

            graftings_ = nullptr;
            attached_ = false;

        }

        virtual ~ForwardingSwarm() {
            delete graftings_;
        }

        ForwardingSwarm(const ForwardingSwarm& other) : // copy constructor
                swarms_(other.swarms_), start_(other.start_), size_(other.size_), mass_(other.mass_), num_different_(other.num_different_),
                num_singletons_(other.num_singletons_), max_rad_(other.max_rad_), graftings_(other.graftings_->clone()),
                attached_(other.attached_) {

            // nothing else to do

        }

        ForwardingSwarm(ForwardingSwarm&& other) noexcept : // move constructor
                swarms_(other.swarms_), start_(other.start_),
                size_(other.size_), mass_(other.mass_), num_different_(other.num_different_), num_singletons_(other.num_singletons_),
                max_rad_(other.max_rad_), graftings_(other.graftings_), attached_(other.attached_) {

            other.graftings_ = nullptr;

        }

        ForwardingSwarm& operator=(const ForwardingSwarm& other) { // copy assignment operator

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

            graftings_ = other.graftings_->clone();
            attached_ = other.attached_;

            return *this;

        }

        ForwardingSwarm& operator=(ForwardingSwarm&& other) noexcept { // move assignment operator

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

        ForwardingSwarm* clone() const override {
            return new ForwardingSwarm(*this);
        }


        inline numSeqs_t seed() const override {
            return swarms_->get_member(start_);
        }
        inline SwarmEntry* seed_entry() const override {
            return new ForwardingSwarmEntry<ForwardingSwarm>(this, 0);
        }
        inline SwarmEntry* seed_entry(SwarmEntry* entry) const override {

            dynamic_cast<ForwardingSwarmEntry<ForwardingSwarm>*>(entry)->reassign(this, 0);
            return entry;

        }

        inline numSeqs_t get(numSeqs_t i) override {
            return swarms_->get_member(start_ + i);
        }
        inline SwarmEntry* get_entry(numSeqs_t i) override {
            return new ForwardingSwarmEntry<ForwardingSwarm>(this, i);
        }
        inline SwarmEntry* get_entry(numSeqs_t i, SwarmEntry* entry) const override {

            dynamic_cast<ForwardingSwarmEntry<ForwardingSwarm>*>(entry)->reassign(this, i);
            return entry;

        }

        inline numSeqs_t member(numSeqs_t i) const override {
            return swarms_->get_member(start_ + i);
        }

        inline numSeqs_t parent(numSeqs_t i) const override {
            return swarms_->get_parent(start_ + i);
        }

        inline numSeqs_t gen(numSeqs_t i) const override {
            return swarms_->get_gen(start_ + i);
        }

        inline dist_t rad(numSeqs_t i) const override {
            return swarms_->get_rad(start_ + i);
        }

        inline dist_t parent_dist(numSeqs_t i) const override {
            return swarms_->get_parent_dist(start_ + i);
        }

        inline numSeqs_t size() const override {
            return size_;
        }

        inline numSeqs_t num_different() const override {
            return num_different_;
        }

        inline numSeqs_t num_singletons() const override {
            return num_singletons_;
        }

        inline numSeqs_t mass() const override {
            return mass_;
        }

        void sort_next_generation(const AmpliconCollection& ac, numSeqs_t p) override {

            std::vector<numSeqs_t> permutation(size() - p);
            std::iota(permutation.begin(), permutation.end(), p);

            std::sort(permutation.begin(), permutation.end(),
                      [&](const numSeqs_t a, const numSeqs_t b) {
                          return (ac.ab(swarms_->get_member(start_ + a)) > ac.ab(swarms_->get_member(start_ + b))) ||
                                 ((ac.ab(swarms_->get_member(start_ + a)) == ac.ab(swarms_->get_member(start_ + b))) &&
                                  (strcmp(ac.id(swarms_->get_member(start_ + a)), ac.id(swarms_->get_member(start_ + b))) < 0));
                      }
            );

            std::vector<bool> done(size());
            for (numSeqs_t i = p; i < size(); ++i) {

                if (!done[i]) {

                    done[i] = true;
                    numSeqs_t prev_j = i;
                    numSeqs_t j = permutation[i - p];

                    while (i != j) {

                        swarms_->swap(start_ + prev_j, start_ + j);
                        done[j] = true;
                        prev_j = j;
                        j = permutation[j - p];

                    }
                }

            }

        }

        inline void append(numSeqs_t m, numSeqs_t g, dist_t r, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) override {

            swarms_->append(m, p, dist, g, r);

            mass_ += ab;
            num_different_ += new_seq;
            num_singletons_ += (ab == 1);

            max_rad_ = std::max(max_rad_, r);

            size_++;

        }

        inline void append(numSeqs_t m, numSeqs_t p, dist_t dist, numSeqs_t ab, bool new_seq) override {

            swarms_->append(m, swarms_->get_member(start_ + p), dist, swarms_->get_gen(start_ + p) + 1, swarms_->get_rad(start_ + p) + dist);

            mass_ += ab;
            num_different_ += new_seq;
            num_singletons_ += (ab == 1);

            max_rad_ = std::max(max_rad_, swarms_->get_rad(start_ + p) + dist);

            size_++;

        }

        void attach(const numSeqs_t p, const numSeqs_t c, Swarm* cs, const dist_t d) override {

            if (graftings_ == nullptr) {
                graftings_ = new SimpleGraftingInfo();
            }

            graftings_->add_link(p, c, cs, d);
            cs->mark_as_attached();

        }

        void mark_as_attached() override {
            attached_ = true;
        }

        bool is_attached() const override {
            return attached_;
        }

        numSeqs_t total_size() const override {

            numSeqs_t sum = size();
            if (graftings_ != nullptr) {

                for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                    sum += graftings_->get_child_swarm(i)->size();
                }

            }

            return sum;

        }

        numSeqs_t total_num_different() const override {

            numSeqs_t sum = num_different();
            if (graftings_ != nullptr) {

                for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                    sum += graftings_->get_child_swarm(i)->num_different();
                }

            }

            return sum;

        }

        numSeqs_t total_num_singletons() const override {

            numSeqs_t sum = num_singletons();
            if (graftings_ != nullptr) {

                for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                    sum += graftings_->get_child_swarm(i)->num_singletons();
                }

            }

            return sum;

        }

        numSeqs_t total_mass() const override {

            numSeqs_t sum = mass();
            if (graftings_ != nullptr) {

                for (numSeqs_t i = 0; i < graftings_->num_links(); i++) {
                    sum += graftings_->get_child_swarm(i)->mass();
                }

            }

            return sum;

        }

        dist_t max_rad() const override {
            return max_rad_;
        }

        lenSeqs_t max_gen() const override {
            return (size_ > 0) ? swarms_->get_gen(start_ + size_ - 1) : 0;
        }

        GraftingInfo* get_grafting_info() const override {
            return graftings_;
        }

    protected:
        ForwardingSwarm() = default;

        // the information on each swarm member are stored separately but in the same order in multiple containers
        SharingSwarms* swarms_; // pool-internal integer id of swarm member
        numSeqs_t start_;

        // overall characteristics of the swarm (not considering the grafted amplicons)
        numSeqs_t size_; // number of members of the swarm
        numSeqs_t mass_; // sum of the abundances of amplicons
        numSeqs_t num_different_; // number of different sequences
        numSeqs_t num_singletons_; // number of singletons
        dist_t max_rad_; // maximum radius

        GraftingInfo* graftings_; // grafting links with this swarm as the parent (nullptr when there are none)
        bool attached_; // Boolean flag indicating whether this swarm is attached (as a child) to another swarm

    };

}

#endif //GEFAST_SWARMS_HPP
