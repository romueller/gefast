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

    protected:
        std::vector<S> swarms_per_pool_; // collection of swarms / clusters per pool

    };

}

#endif //GEFAST_SWARMSTORAGES_HPP
