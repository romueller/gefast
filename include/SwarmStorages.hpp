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

#ifndef GEFAST_SWARMSTORAGES_HPP
#define GEFAST_SWARMSTORAGES_HPP

#include "Base.hpp"
#include "Swarms.hpp"

namespace GeFaST {

    /*
     * Class for managing swarm collections coming from one or more AmpliconCollection instances (pools).
     *
     * Basic implementation storing the swarms / clusters as separate entities.
     */
    template<typename S>
    class SimpleSwarmStorage : public SwarmStorage {

    public:
        explicit SimpleSwarmStorage(const numSeqs_t num_pools) {
            swarms_per_pool_ = std::vector<S>(num_pools);
        }

        virtual ~SimpleSwarmStorage() = default;

        SimpleSwarmStorage(const SimpleSwarmStorage& other) = default; // copy constructor

        SimpleSwarmStorage(SimpleSwarmStorage&& other) = default; // move constructor

        SimpleSwarmStorage& operator=(const SimpleSwarmStorage& other) = default; // copy assignment operator

        SimpleSwarmStorage& operator=(SimpleSwarmStorage&& other) = default; // move assignment operator

        SimpleSwarmStorage* clone() const override { // deep-copy clone method
            return new SimpleSwarmStorage(*this);
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

            return sizeof(SimpleSwarmStorage) // SimpleSwarmStorage itself
                + allocated_swarms // swarms_per_pool_ (actually existing swarm collections)
                + sizeof(S) * (swarms_per_pool_.capacity() - swarms_per_pool_.size()); // swarms_per_pool_ (allocated, empty swarm collections in vector)

        }

        void show_memory() const override {

            size_t allocated_swarms = 0;
            for (auto& s : swarms_per_pool_) {
                allocated_swarms += s.size_in_bytes();
            }

            std::cout << "##################################################" << std::endl;
            std::cout << "# SimpleSwarmStorage" << std::endl;
            std::cout << "# " << std::endl;
            std::cout << "# sizeof(SimpleSwarmStorage<S>): " << sizeof(SimpleSwarmStorage<S>) << " bytes" << std::endl;
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

}

#endif //GEFAST_SWARMSTORAGES_HPP
