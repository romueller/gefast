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

#ifndef GEFAST_AMPLICONSTORAGES_HPP
#define GEFAST_AMPLICONSTORAGES_HPP

#include "Base.hpp"
#include "space/Basics.hpp"

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

        size_t size_in_bytes() const override {

            size_t allocated_pools = 0;
            for (auto& p : pools_) {
                allocated_pools += p.size_in_bytes();
            }

            return sizeof(LengthPoolsAmpliconStorage) // LengthPoolsAmpliconStorage itself
                + allocated_pools // pools_ (actually existing amplicon collections)
                + sizeof(C) * (pools_.capacity() - pools_.size()) // pools_ (allocated, empty amplicon collections in vector)
                + map_pod_size_in_bytes(pool_map_) - sizeof(std::map<lenSeqs_t, numSeqs_t>); // pool_map_ (do not count size of map object twice)

        }

        void show_memory() const override {

            size_t allocated_pools = 0;
            for (auto& p : pools_) {
                allocated_pools += p.size_in_bytes();
            }

            std::cout << "##################################################" << std::endl;
            std::cout << "# LengthPoolsAmpliconStorage" << std::endl;
            std::cout << "#  - number of pools: " << pools_.size() << std::endl;
            std::cout << "#  - max_len_: " << max_len_ << std::endl;
            std::cout << "# " << std::endl;
            std::cout << "# sizeof(LengthPoolsAmpliconStorage): " << sizeof(LengthPoolsAmpliconStorage) << " bytes" << std::endl;
            std::cout << "# sizeof(std::vector<C>): " << sizeof(std::vector<C>) << " bytes" << std::endl;
            std::cout << "# sizeof(lenSeqs_t): " << sizeof(lenSeqs_t) << " bytes" << std::endl;
            std::cout << "# sizeof(std::map<lenSeqs_t, numSeqs_t>): " << sizeof(std::map<lenSeqs_t, numSeqs_t>) << " bytes" << std::endl;
            std::cout << "# " << std::endl;
            std::cout << "# Total size: " << size_in_bytes() << " bytes" << std::endl;
            std::cout << "# pools_ (actual amplicon collections): " << allocated_pools << " bytes" << std::endl;
            std::cout << "# pools_ (additional reserved space in vector): " << sizeof(C) * (pools_.capacity() - pools_.size()) << " bytes" << std::endl;
            std::cout << "# pool_map_: " << map_pod_size_in_bytes(pool_map_) << " bytes" << std::endl;
            std::cout << "# " << std::endl;
            std::cout << "# Size of amplicon collections:" << std::endl;
            for (auto p = 0; p < num_pools(); p++) std::cout << "# " << p << ": " << get_pool(p).size_in_bytes() << " bytes" << std::endl;
            std::cout << "##################################################" << std::endl;

        }

    protected:
        LengthPoolsAmpliconStorage() = default;

        std::vector<C> pools_; // collection of amplicon pools
        lenSeqs_t max_len_ = 0; // maximum length of the amplicons in the pools

        std::map<lenSeqs_t, numSeqs_t> pool_map_; // mapping between sequence length and pools

    };


    /*
     * Class for managing amplicons stored in a single AmpliconCollection instance (pool)
     * representing overall collection of amplicons.
     *
     * Instances are constructed via an AmpliconStorageFactory.
     */
    template<typename C>
    class SimpleFeatureAmpliconStorage : public FeatureAmpliconStorage {

        friend class AmpliconStorageFactory;

    public:
        SimpleFeatureAmpliconStorage(const C& pool) : pool_(pool) {
            // nothing else to do
        }

        SimpleFeatureAmpliconStorage* clone() const override { // deep-copy clone method
            return new SimpleFeatureAmpliconStorage(*this);
        }


        void add_amplicon(const Defline& dl, const std::string& seq, const ExtraInfo& extra) override {

            pool_.emplace_back(dl, seq, extra);
            max_len_ = std::max(max_len_, seq.length());

        }

        FeatureAmpliconCollection& get_pool(const numSeqs_t i) override {
            return pool_;
        }

        const FeatureAmpliconCollection& get_pool(const numSeqs_t i) const override {
            return pool_;
        }

        numSeqs_t num_pools() const override {
            return 1;
        }

        numSeqs_t num_amplicons() const override {
            return pool_.size();
        }

        numSeqs_t size_pool(numSeqs_t i) const override {
            return pool_.size();
        }

        lenSeqs_t max_length() const override {
            return max_len_;
        }

        void finalise() override {
            pool_.sort();
        }

        void print(const std::string& file, const Configuration& config) const override {

            // clear file if it already exists
            std::ofstream out_stream(file);
            out_stream.close();

            // print contents of the pool in one file
            pool_.print(file, config);

        }

        size_t size_in_bytes() const override {

            std::cerr << "ERROR: SimpleFeatureAmpliconStorage is currently not part of the space-efficiency analysis." << std::endl;
            return 0;

        }

        void show_memory() const override {
            std::cerr << "ERROR: SimpleFeatureAmpliconStorage is currently not part of the space-efficiency analysis." << std::endl;
        }

    protected:
        SimpleFeatureAmpliconStorage() = default;

        C pool_; // single amplicon underlying the amplicon storage
        lenSeqs_t max_len_ = 0; // maximum length of the amplicons in the pool

    };

}

#endif //GEFAST_AMPLICONSTORAGES_HPP
