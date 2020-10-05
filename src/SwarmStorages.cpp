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

#include "../include/SwarmStorages.hpp"

namespace GeFaST {

    /* === SimpleSwarmStorage === */

    SimpleSwarmStorage::SimpleSwarmStorage(const numSeqs_t num_pools) {
        swarms_per_pool_ = std::vector<SimpleSwarms<SimpleSwarm>>(num_pools);
    }

    SimpleSwarmStorage* SimpleSwarmStorage::clone() const {
        return new SimpleSwarmStorage(*this);
    }
    
    Swarms& SimpleSwarmStorage::get_swarms(numSeqs_t p) {
        return swarms_per_pool_[p];
    }

    const Swarms& SimpleSwarmStorage::get_swarms(numSeqs_t p) const {
        return swarms_per_pool_[p];
    }

    numSeqs_t SimpleSwarmStorage::num_pools() const {
        return swarms_per_pool_.size();
    }

    numSeqs_t SimpleSwarmStorage::num_swarms() const {

        numSeqs_t sum = 0;
        for (auto& sws : swarms_per_pool_) {
            sum += sws.size();
            for (auto i = 0; i < sws.size(); i++) {
                sum -= sws.get(i).is_attached();
            }
        }

        return sum;

    }

    numSeqs_t SimpleSwarmStorage::num_swarms(numSeqs_t p) const {
        return swarms_per_pool_[p].size();
    }

}