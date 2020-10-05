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

#ifndef GEFAST_AMPLICONSTORAGES_HPP
#define GEFAST_AMPLICONSTORAGES_HPP

#include "Base.hpp"

namespace GeFaST {

    /*
     * Class for managing amplicons stored one or more AmpliconCollection instances (pools),
     * which are separated from each other by sufficiently large gaps in the length distribution
     * of the overall collection of amplicons.
     * These gaps ensure that distances considering the number of edit operations cannot find
     * similar amplicons across pools.
     *
     * Instances are constructed via an AmpliconStorageFactory.
     */
    template<typename C>
    class LengthPoolsAmpliconStorage : public AmpliconStorage {

        friend class AmpliconStorageFactory;

    public:
        LengthPoolsAmpliconStorage* clone() const override { // deep-copy clone method
            return new LengthPoolsAmpliconStorage(*this);
        }


        void add_amplicon(const Defline& dl, const std::string& seq, const ExtraInfo& extra) override {

            pools_[pool_map_[seq.length()]].emplace_back(dl, seq, extra);
            max_len_ = std::max(max_len_, seq.length());

        }

        AmpliconCollection& get_pool(const numSeqs_t i) override {
            return pools_[i];
        }

        const AmpliconCollection& get_pool(const numSeqs_t i) const override {
            return pools_[i];
        }

        numSeqs_t num_pools() const override {
            return pools_.size();
        }

        numSeqs_t num_amplicons() const override {

            numSeqs_t num = 0;
            for (auto& p : pools_) {
                num += p.size();
            }

            return num;

        }

        numSeqs_t size_pool(numSeqs_t i) const override {
            return pools_[i].size();
        }

        lenSeqs_t max_length() const override {
            return max_len_;
        }

        void finalise() override {

            for (numSeqs_t i = 0; i < pools_.size(); i++) {
                pools_[i].sort();
            }

            std::map<lenSeqs_t, numSeqs_t>().swap(pool_map_);

        }

        void print(const std::string& file, const Configuration& config) const override {

            // clear file if it already exists
            std::ofstream out_stream(file);
            out_stream.close();

            // print contents of all pools in one file
            for (numSeqs_t i = 0; i < pools_.size(); i++) {
                pools_[i].print(file, config);
            }

        }

    protected:
        LengthPoolsAmpliconStorage() = default;

        std::vector<C> pools_; // collection of amplicon pools
        lenSeqs_t max_len_ = 0; // maximum length of the amplicons in the pools

        std::map<lenSeqs_t, numSeqs_t> pool_map_; // mapping between sequence length and pools

    };

}

#endif //GEFAST_AMPLICONSTORAGES_HPP
