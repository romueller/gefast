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

#ifndef GEFAST_SWARMREPRESENTATIONS_HPP
#define GEFAST_SWARMREPRESENTATIONS_HPP

#include "../Swarms.hpp"
#include "NumericArrays.hpp"

namespace GeFaST {

    /*
     * Class for managing swarm collections coming from one or more AmpliconCollection instances (pools).
     *
     * Basic implementation storing the swarms / clusters as separate entities.
     */
    template<typename S>
    class DacsSwarmStorage : public SwarmStorage {

    public:
        explicit DacsSwarmStorage(const AmpliconStorage& as, const SpaceLevenshteinConfiguration& config) {
            for (auto i = 0; i < as.num_pools(); i++) {
                swarms_per_pool_.emplace_back(as.size_pool(i), config.chunk_lengths);
            }

        }

        virtual ~DacsSwarmStorage() = default;

        DacsSwarmStorage(const DacsSwarmStorage& other) = default; // copy constructor

        DacsSwarmStorage(DacsSwarmStorage&& other) noexcept = default; // move constructor

        DacsSwarmStorage& operator=(const DacsSwarmStorage& other) = default; // copy assignment operator

        DacsSwarmStorage& operator=(DacsSwarmStorage&& other) noexcept = default; // move assignment operator

        DacsSwarmStorage* clone() const override { // deep-copy clone method
            return new DacsSwarmStorage(*this);
        }


        Swarms& get_swarms(numSeqs_t p) override {
            return swarms_per_pool_[p];
        }

        const Swarms& get_swarms(numSeqs_t p) const override {
            return swarms_per_pool_[p];
        }

        numSeqs_t num_pools() const override {
            return swarms_per_pool_.size();
        }

        numSeqs_t num_swarms() const override {

            numSeqs_t sum = 0;
            for (auto& sws : swarms_per_pool_) {
                sum += sws.size();
                for (auto i = 0; i < sws.size(); i++) {
                    sum -= sws.get(i).is_attached();
                }
            }

            return sum;

        }

        numSeqs_t num_swarms(numSeqs_t p) const override {
            return swarms_per_pool_[p].size();
        }

        size_t size_in_bytes() const override {

            size_t allocated_swarms = 0;
            for (auto& s : swarms_per_pool_) {
                allocated_swarms += s.size_in_bytes();
            }

            return sizeof(DacsSwarmStorage) // DacsSwarmStorage itself
                + allocated_swarms // swarms_per_pool_ (actually existing swarm collections)
                + sizeof(S) * (swarms_per_pool_.capacity() - swarms_per_pool_.size()); // swarms_per_pool_ (allocated, empty swarm collections in vector)

        }

        void show_memory() const override {

            size_t allocated_swarms = 0;
            for (auto& s : swarms_per_pool_) {
                allocated_swarms += s.size_in_bytes();
            }

            std::cout << "##################################################" << std::endl;
            std::cout << "# DacsSwarmStorage" << std::endl;
            std::cout << "# " << std::endl;
            std::cout << "# sizeof(DacsSwarmStorage<S>): " << sizeof(DacsSwarmStorage<S>) << " bytes" << std::endl;
            std::cout << "# sizeof(std::vector<S>): " << sizeof(std::vector<S>) << " bytes" << std::endl;
            std::cout << "# sizeof(S): " << sizeof(S) << " bytes" << std::endl;
            std::cout << "# " << std::endl;
            std::cout << "# Total size: " << size_in_bytes() << " bytes" << std::endl;
            std::cout << "# swarms_per_pool_ (actual swarm collections): " << allocated_swarms << " bytes" << std::endl;
            std::cout << "# swarms_per_pool_ (additional reserved space in vector): " << sizeof(S) * (swarms_per_pool_.capacity() - swarms_per_pool_.size()) << " bytes" << std::endl;
            std::cout << "# " << std::endl;
            std::cout << "# Size of swarm collections:" << std::endl;
            for (auto p = 0; p < num_pools(); p++) std::cout << "# " << p << ": " << get_swarms(p).size_in_bytes() << " bytes" << std::endl;
            std::cout << "##################################################" << std::endl;

        }

    protected:
        std::vector<S> swarms_per_pool_; // collection of swarms / clusters per pool

    };


    /*
     * The information of each swarm member are split and stored separately (but in the same order)
     * in different DACs sequences.
     */
    class SharingSplitDacsSwarms : public Swarms {

    public:
        SharingSplitDacsSwarms() = default;

        SharingSplitDacsSwarms(const numSeqs_t pool_size, const ChunkLengths& config);

        Swarm& initialise_cluster(numSeqs_t s, numSeqs_t ab) override;

        ~SharingSplitDacsSwarms();

        SharingSplitDacsSwarms(const SharingSplitDacsSwarms& other); // copy constructor

        SharingSplitDacsSwarms(SharingSplitDacsSwarms&& other) noexcept; // move constructor

        SharingSplitDacsSwarms& operator=(const SharingSplitDacsSwarms& other); // copy assignment operator

        SharingSplitDacsSwarms& operator=(SharingSplitDacsSwarms&& other) noexcept; // move assignment operator

        SharingSplitDacsSwarms* clone() const override;


        numSeqs_t size() const override;

        Swarm& get(const numSeqs_t i) override;

        const Swarm& get(const numSeqs_t i) const override;

        void clear() override;

        void finalise() override;

        void append(numSeqs_t m, numSeqs_t p, dist_t d, numSeqs_t g, dist_t r);

        void swap(numSeqs_t i, numSeqs_t j);

        void sort_next_generation(const AmpliconCollection& ac);

        size_t size_in_bytes() const override;

        void show_memory(numSeqs_t pid) const override;


    private:
        typedef SplitDacs<numSeqs_t, size_t> dacs_t;
        typedef SplitDacsDirectBuilder<numSeqs_t, size_t> builder_t;

        std::vector<Swarm*> swarms_; // swarms contained in the collection

        // the information on each swarm member,
        // during clustering except for members of current swarm in current and next generation
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> all_members_; // pool-internal integer ids of swarm members
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> all_parents_; // pool-internal integer ids of its parents
        DacsCapacityArray<dacs_t, builder_t, lenSeqs_t> all_parent_dists_; // distances to parent
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> all_gens_; // generation numbers of members
        DacsCapacityArray<dacs_t, builder_t, lenSeqs_t> all_rads_; // radii of members

        // information on members of current generation
        std::vector<numSeqs_t> cur_members_, cur_parents_; // pool-internal integer ids of members and parents
        std::vector<lenSeqs_t> cur_dists_, cur_rads_; // distances to parent and radii
        numSeqs_t cur_gen_; // generation number

        // information on members of next generation
        std::vector<numSeqs_t> next_members_, next_parents_; // pool-internal integer ids of members and parents
        std::vector<lenSeqs_t> next_dists_, next_rads_; // distances to parent and radii

        numSeqs_t total_num_members_; // current number of entries in the DACs sequences all_members_ etc.
        bool finalised_; // flag indicating whether the construction process has already been finished


        /*
         * Helper method for finishing the work on the current generation of subseeds.
         * Moves the information on the members of the current generation to the DACs sequences.
         * Prepares cur_members_ etc. and next_members_ etc. to extend the cluster further.
         */
        void transfer_generations();

        /*
         * Helper methods used by below Swarm representation to access the necessary member information during
         * construction and later on.
         */
        numSeqs_t member(numSeqs_t i) const;
        numSeqs_t parent(numSeqs_t i) const;
        numSeqs_t gen(numSeqs_t i) const;
        dist_t rad(numSeqs_t i) const;
        dist_t parent_dist(numSeqs_t i) const;



        /*
         * Representation of a swarm storing only the overall characteristics (e.g. mass).
         * The information on each member of the swarm (parent, generation etc.) are stored in SharingSplitDacsSwarms instance.
         * Requests of such information are forwarded to the SharingSplitDacsSwarms instance.
         * Swarms grafted to a particular swarm instance are only added implicitly
         * through information on the grafting links.
         */
        class ForwardingDacsSwarm : public Swarm {
        public:

            /*
             * Create new swarm and initialise it with member m of abundance ab as the seed.
             */
            ForwardingDacsSwarm(SharingSplitDacsSwarms& as, numSeqs_t ab);

            virtual ~ForwardingDacsSwarm();

            ForwardingDacsSwarm(const ForwardingDacsSwarm& other); // copy constructor

            ForwardingDacsSwarm(ForwardingDacsSwarm&& other) noexcept; // move constructor

            ForwardingDacsSwarm& operator=(const ForwardingDacsSwarm& other); // copy assignment operator

            ForwardingDacsSwarm& operator=(ForwardingDacsSwarm&& other) noexcept; // move assignment operator

            ForwardingDacsSwarm* clone() const override;


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

            size_t size_in_bytes() const override;

        protected:
            ForwardingDacsSwarm() = default;

            // the information on each swarm member are stored separately but in the same order in multiple containers
            SharingSplitDacsSwarms* swarms_; // Swarms representation containing all member information
            numSeqs_t start_; // start position of this Swarm in the collective storage

            // overall characteristics of the swarm (not considering the grafted amplicons)
            numSeqs_t size_; // number of members of the swarm
            numSeqs_t mass_; // sum of the abundances of amplicons
            numSeqs_t num_different_; // number of different sequences
            numSeqs_t num_singletons_; // number of singletons
            dist_t max_rad_; // maximum radius

            GraftingInfo* graftings_; // grafting links with this swarm as the parent (nullptr when there are none)
            bool attached_; // Boolean flag indicating whether this swarm is attached (as a child) to another swarm

        };

    };


    /*
     * The information of each swarm member are split and stored separately (but in the same order)
     * in different DACs sequences.
     * In addition, the (numeric) overall characteristics of the swarms are also stored as DACs sequences.
     */
    class SharingSplitFullDacsSwarms : public Swarms {

    public:
        SharingSplitFullDacsSwarms() = default;

        SharingSplitFullDacsSwarms(const numSeqs_t pool_size, const ChunkLengths& config);

        Swarm& initialise_cluster(numSeqs_t s, numSeqs_t ab) override;

        ~SharingSplitFullDacsSwarms();

        SharingSplitFullDacsSwarms(const SharingSplitFullDacsSwarms& other); // copy constructor

        SharingSplitFullDacsSwarms(SharingSplitFullDacsSwarms&& other) noexcept; // move constructor

        SharingSplitFullDacsSwarms& operator=(const SharingSplitFullDacsSwarms& other); // copy assignment operator

        SharingSplitFullDacsSwarms& operator=(SharingSplitFullDacsSwarms&& other) noexcept; // move assignment operator

        SharingSplitFullDacsSwarms* clone() const override;


        numSeqs_t size() const override;

        Swarm& get(const numSeqs_t i) override;

        const Swarm& get(const numSeqs_t i) const override;

        void clear() override;

        void finalise() override;

        void append(numSeqs_t m, numSeqs_t p, dist_t d, numSeqs_t g, dist_t r, numSeqs_t ab, bool new_seq);

        void swap(numSeqs_t i, numSeqs_t j);

        void sort_next_generation(const AmpliconCollection& ac);

        size_t size_in_bytes() const override;

        void show_memory(numSeqs_t pid) const override;


    private:
        typedef SplitDacs<numSeqs_t, size_t> dacs_t;
        typedef SplitDacsDirectBuilder<numSeqs_t, size_t> builder_t;

        std::vector<Swarm*> swarms_; // swarms contained in the collection

        // the information on each swarm member,
        // during clustering except for members of current swarm in current and next generation
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> all_members_; // pool-internal integer ids of swarm members
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> all_parents_; // pool-internal integer ids of its parents
        DacsCapacityArray<dacs_t, builder_t, lenSeqs_t> all_parent_dists_; // distances to parent
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> all_gens_; // generation numbers of members
        DacsCapacityArray<dacs_t, builder_t, lenSeqs_t> all_rads_; // radii of members

        // information on members of current generation
        std::vector<numSeqs_t> cur_members_, cur_parents_; // pool-internal integer ids of members and parents
        std::vector<lenSeqs_t> cur_dists_, cur_rads_; // distances to parent and radii
        numSeqs_t cur_gen_; // generation number

        // information on members of next generation
        std::vector<numSeqs_t> next_members_, next_parents_; // pool-internal integer ids of members and parents
        std::vector<lenSeqs_t> next_dists_, next_rads_; // distances to parent and radii

        // information on the swarms,
        // during construction except for the current swarm
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> swarms_start_; // start position of each swarm in all_members_ etc.
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> swarms_size_; // number of members of each swarm
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> swarms_mass_; // sum of the abundances of amplicons in each swarm
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> swarms_num_different_; // number of different sequences in each swarm
        DacsCapacityArray<dacs_t, builder_t, numSeqs_t> swarms_num_singletons_; // number of singletons in each swarm
        DacsCapacityArray<dacs_t, builder_t, lenSeqs_t> swarms_max_rad_; // maximum radius of each swarm
        std::vector<GraftingInfo*> swarms_graftings_; // grafting links of each swarm
        std::vector<bool> swarms_attached_; // flags indicating whether the swarms are attached

        // values of currently built swarm
        numSeqs_t cur_swarm_start_; // start position of the swarm
        numSeqs_t cur_swarm_size_; // number of members of the swarm
        numSeqs_t cur_swarm_mass_; // sum of the abundances of amplicons
        numSeqs_t cur_swarm_num_different_; // number of different sequences
        numSeqs_t cur_swarm_num_singletons_; // number of singletons
        dist_t cur_swarm_max_rad_; // maximum radius

        numSeqs_t total_num_members_; // current number of entries in the DACs sequences all_members_ etc.
        bool finalised_; // flag indicating whether the construction process has already been finished


        /*
         * Helper method for finishing the work on the current generation of subseeds.
         * Moves the information on the members of the current generation to the DACs sequences.
         * Prepares cur_members_ etc. and next_members_ etc. to extend the cluster further.
         */
        void transfer_generations();

        /*
         * Helper method for finishing the work on the current swarm.
         * Moves the information on the swarm to the DACs sequences.
         * Prepares cur_swarm_start_ etc. for the next swarm.
         */
        void close_swarm(numSeqs_t s, numSeqs_t ab);

        /*
         * Helper methods used by below Swarm representation to access the necessary member and swarm information during
         * construction and later on.
         */
        numSeqs_t member(numSeqs_t sid, numSeqs_t i) const;
        numSeqs_t parent(numSeqs_t sid, numSeqs_t i) const;
        numSeqs_t gen(numSeqs_t sid, numSeqs_t i) const;
        dist_t rad(numSeqs_t sid, numSeqs_t i) const;
        dist_t parent_dist(numSeqs_t sid, numSeqs_t i) const;

        numSeqs_t swarm_size(numSeqs_t sid) const;
        numSeqs_t swarm_num_different(numSeqs_t sid) const;
        numSeqs_t swarm_num_singletons(numSeqs_t sid) const;
        numSeqs_t swarm_mass(numSeqs_t sid) const;



        /*
         * Representation of a swarm storing only an internal swarm id.
         * The information on each member of the swarm (parent, generation etc.) as well as the overall characteristics
         * of the swarm are stored in the SharingSplitFullDacsSwarms instance.
         * Requests of such information are forwarded to the SharingSwarms instance.
         * Swarms grafted to a particular swarm instance are only added implicitly
         * through information on the grafting links.
         */
        class ForwardingDacsSwarm : public Swarm {
        public:

            /*
             * Create new swarm and initialise it with member m of abundance ab as the seed.
             */
            ForwardingDacsSwarm(SharingSplitFullDacsSwarms* swarms, numSeqs_t ab);

            virtual ~ForwardingDacsSwarm() = default;

            ForwardingDacsSwarm(const ForwardingDacsSwarm& other) = default; // copy constructor

            ForwardingDacsSwarm(ForwardingDacsSwarm&& other) noexcept = default; // move constructor

            ForwardingDacsSwarm& operator=(const ForwardingDacsSwarm& other) = default; // copy assignment operator

            ForwardingDacsSwarm& operator=(ForwardingDacsSwarm&& other) noexcept = default; // move assignment operator

            ForwardingDacsSwarm* clone() const override;


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

            size_t size_in_bytes() const override;

        protected:
            ForwardingDacsSwarm() = default;

            SharingSplitFullDacsSwarms* swarms_; // Swarms representation containing all information
            numSeqs_t sid_; // internal id of this swarm in swarms_

        };

    };

}

#endif //GEFAST_SWARMREPRESENTATIONS_HPP
