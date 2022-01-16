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

#include <iostream>

#include "../include/Factories.hpp"
#include "../include/modes/DereplicationMode.hpp"
#include "../include/AmpliconCollections.hpp"
#include "../include/AmpliconStorages.hpp"
#include "../include/Distances.hpp"
#include "../include/FeatureBuilders.hpp"
#include "../include/Preprocessors.hpp"
#include "../include/Clusterers.hpp"
#include "../include/ClusterRefiners.hpp"
#include "../include/OutputGenerators.hpp"
#include "../include/AuxiliaryData.hpp"
#include "../include/SwarmStorages.hpp"
#include "../include/modes/PreprocessingOnlyMode.hpp"
#include "../include/modes/QualityAlignmentScoreMode.hpp"
#include "../include/modes/ConsistencyMode.hpp"
#include "../include/QualityWeightedDistances.hpp"
#include "../include/SpacePartitioning.hpp"
#include "../include/modes/SpaceLevenshteinMode.hpp"
#include "../include/space/FlexibleAmpliconCollection.hpp"
#include "../include/space/FlexibleRelations.hpp"
#include "../include/space/NumericArrays.hpp"
#include "../include/space/StringArrays.hpp"
#include "../include/space/SwarmRepresentations.hpp"
#include "../include/space/k2trees/BasicK2Tree.hpp"
#include "../include/space/k2trees/HybridK2Tree.hpp"
#include "../include/space/k2trees/NaiveK2Tree.hpp"
#include "../include/space/k2trees/PartitionedRectangularK2Tree.hpp"
#include "../include/space/k2trees/PartitionedVariableK2Tree.hpp"
#include "../include/space/k2trees/RectangularK2Tree.hpp"

namespace GeFaST {

    Configuration* ConfigurationFactory::create(int argc, const char* argv[]) {

        std::string mode = argv[1];
        Configuration* config = nullptr;

        switch (get_configuration_option(mode)) {
            case CF_PREPROCESSING_ONLY:
                config = new PreprocessingOnlyConfiguration(argc, argv);
                break;

            case CF_DEREPLICATION:
                config = new DereplicationConfiguration(argc, argv);
                break;

            case CF_LEVENSHTEIN:
                config = new LevenshteinConfiguration(argc, argv);
                break;

            case CF_QUALITY_LEVENSHTEIN:
                config = new QualityLevenshteinConfiguration(argc, argv);
                break;

            case CF_ALIGNMENT_SCORE:
                config = new AlignmentScoreConfiguration(argc, argv);
                break;

            case CF_QUALITY_ALIGNMENT_SCORE:
                config = new QualityAlignmentScoreConfiguration(argc, argv);
                break;

            case CF_ALIGNMENT_FREE:
                config = new AlignmentFreeConfiguration(argc, argv);
                break;

            case CF_CONSISTENCY:
                config = new ConsistencyConfiguration(argc, argv);
                break;

            case CF_SPACE_LEVENSHTEIN:
                config = new SpaceLevenshteinConfiguration(argc, argv);
                break;

            case CF_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown mode '" << mode << "'." << std::endl;

        }

        return config;

    }



    AmpliconStorage* AmpliconStorageFactory::create(const AmpliconStorageOption opt, const DataStatistics<>& ds,
            const Configuration& config) {

        AmpliconStorage* res = nullptr;

        switch (opt) {

            // Create multiple pools when there are sufficient gaps in the length distribution,
            // across which amplicons cannot be similar.
            // The pools are not initialised with a particular capacity but grow as the amplicons are inserted.
            case AS_SIMPLE_LENGTH_POOLS: {
                auto tmp = new LengthPoolsAmpliconStorage<SimpleAmpliconCollection>();

                std::vector<lenSeqs_t> lengths = ds.get_all_lengths();
                auto pooling_threshold = static_cast<lenSeqs_t>(std::max(config.main_threshold, config.refinement_threshold));

                lenSeqs_t last_len = ds.get_min_length();
                numSeqs_t pool_id = 0;

                for (auto len : lengths) {

                    if ((last_len + pooling_threshold) < len || len == 0) { // new pool

                        pool_id++;
                        tmp->pools_.emplace_back();

                    }

                    tmp->pool_map_[len] = pool_id;
                    last_len = len;

                }

                tmp->pools_.emplace_back();
                res = tmp;
                break;
            }

            // Create multiple pools when there are sufficient gaps in the length distribution,
            // across which amplicons cannot be similar.
            // Each pool is initialised with a sufficient capacity for the amplicons expected to be inserted.
            case AS_PREPARED_LENGTH_POOLS: {
                auto tmp = new LengthPoolsAmpliconStorage<ArrayAmpliconCollection>();

                std::vector<lenSeqs_t> lengths = ds.get_all_lengths();
                auto pooling_threshold = static_cast<lenSeqs_t>(std::max(config.main_threshold, config.refinement_threshold));

                lenSeqs_t last_len = ds.get_min_length();
                numSeqs_t pool_id = 0;

                numSeqs_t cnt = 0;
                unsigned long long total_len_sequences = 0;
                unsigned long long total_len_headers = 0;

                for (auto len : lengths) {

                    if ((last_len + pooling_threshold) < len) { // new pool

                        pool_id++;
                        tmp->pools_.emplace_back(cnt, total_len_headers, total_len_sequences);

                        cnt = 0;
                        total_len_sequences = 0;
                        total_len_headers = 0;

                    }

                    cnt += ds.get_num_per_length(len);
                    total_len_sequences += (len + 1) * ds.get_num_per_length(len);
                    total_len_headers += ds.get_total_length_headers_per_length(len) + ds.get_num_per_length(len);

                    tmp->pool_map_[len] = pool_id;
                    last_len = len;

                }

                tmp->pools_.emplace_back(cnt, total_len_headers, total_len_sequences);
                res = tmp;
                break;
            }

            // Create multiple pools when there are sufficient gaps in the length distribution,
            // across which amplicons cannot be similar.
            // Each pool is initialised with a sufficient capacity for the amplicons (including their quality scores)
            // expected to be inserted.
            case AS_PREPARED_LENGTH_POOLS_QUALITY: {
                auto tmp = new LengthPoolsAmpliconStorage<ArrayQualityAmpliconCollection>();

                std::vector<lenSeqs_t> lengths = ds.get_all_lengths();
                auto pooling_threshold = static_cast<lenSeqs_t>(std::max(config.main_threshold, config.refinement_threshold));

                lenSeqs_t last_len = ds.get_min_length();
                numSeqs_t pool_id = 0;

                numSeqs_t cnt = 0;
                unsigned long long total_len_sequences = 0;
                unsigned long long total_len_headers = 0;

                for (auto len : lengths) {

                    if ((last_len + pooling_threshold) < len) { // new pool

                        pool_id++;
                        tmp->pools_.emplace_back(cnt, total_len_headers, total_len_sequences);

                        cnt = 0;
                        total_len_sequences = 0;
                        total_len_headers = 0;

                    }

                    cnt += ds.get_num_per_length(len);
                    total_len_sequences += (len + 1) * ds.get_num_per_length(len);
                    total_len_headers += ds.get_total_length_headers_per_length(len) + ds.get_num_per_length(len);

                    tmp->pool_map_[len] = pool_id;
                    last_len = len;

                }

                tmp->pools_.emplace_back(cnt, total_len_headers, total_len_sequences);
                res = tmp;
                break;
            }

            // Create a single pool for all amplicons.
            // The pool is initialised with a sufficient capacity for the amplicons expected to be inserted.
            case AS_PREPARED_POOLS: {
                auto tmp = new LengthPoolsAmpliconStorage<ArrayAmpliconCollection>();

                std::vector<lenSeqs_t> lengths = ds.get_all_lengths();
                numSeqs_t cnt = 0;
                unsigned long long total_len_sequences = 0;
                unsigned long long total_len_headers = 0;

                for (auto len : lengths) {

                    cnt += ds.get_num_per_length(len);
                    total_len_sequences += (len + 1) * ds.get_num_per_length(len);
                    total_len_headers += ds.get_total_length_headers_per_length(len) + ds.get_num_per_length(len);

                    tmp->pool_map_[len] = 0;

                }

                tmp->pools_.emplace_back(cnt, total_len_headers, total_len_sequences);
                res = tmp;
                break;
            }

            // Create a single pool for all amplicons.
            // The pool is initialised with a sufficient capacity for the amplicons (including their quality scores)
            // expected to be inserted.
            case AS_PREPARED_POOLS_QUALITY: {
                auto tmp = new LengthPoolsAmpliconStorage<ArrayQualityAmpliconCollection>();

                std::vector<lenSeqs_t> lengths = ds.get_all_lengths();
                numSeqs_t cnt = 0;
                unsigned long long total_len_sequences = 0;
                unsigned long long total_len_headers = 0;

                for (auto len : lengths) {

                    cnt += ds.get_num_per_length(len);
                    total_len_sequences += (len + 1) * ds.get_num_per_length(len);
                    total_len_headers += ds.get_total_length_headers_per_length(len) + ds.get_num_per_length(len);

                    tmp->pool_map_[len] = 0;

                }

                tmp->pools_.emplace_back(cnt, total_len_headers, total_len_sequences);
                res = tmp;
                break;
            }

            // Create multiple pools when there are sufficient gaps in the length distribution,
            // across which amplicons cannot be similar.
            // Each pool is initialised with a sufficient capacity for the amplicons (including their q-gram profiles)
            // expected to be inserted.
            case AS_PREPARED_QGRAM_LENGTH_POOLS: {
                auto tmp = new LengthPoolsAmpliconStorage<ArrayQgramAmpliconCollection>();

                std::vector<lenSeqs_t> lengths = ds.get_all_lengths();
                auto pooling_threshold = static_cast<lenSeqs_t>(std::max(config.main_threshold, config.refinement_threshold));

                lenSeqs_t last_len = ds.get_min_length();
                numSeqs_t pool_id = 0;

                numSeqs_t cnt = 0;
                unsigned long long total_len_sequences = 0;
                unsigned long long total_len_headers = 0;

                for (auto len : lengths) {

                    if ((last_len + pooling_threshold) < len) { // new pool

                        pool_id++;
                        tmp->pools_.emplace_back(cnt, total_len_headers, total_len_sequences);

                        cnt = 0;
                        total_len_sequences = 0;
                        total_len_headers = 0;

                    }

                    cnt += ds.get_num_per_length(len);
                    total_len_sequences += (len + 1) * ds.get_num_per_length(len);
                    total_len_headers += ds.get_total_length_headers_per_length(len) + ds.get_num_per_length(len);

                    tmp->pool_map_[len] = pool_id;
                    last_len = len;

                }

                tmp->pools_.emplace_back(cnt, total_len_headers, total_len_sequences);
                res = tmp;
                break;
            }

            // Create multiple pools when there are sufficient gaps in the length distribution,
            // across which amplicons cannot be similar.
            // Each pool is initialised with a sufficient capacity for the amplicons (including their quality scores
            // and q-gram profiles) expected to be inserted.
            case AS_PREPARED_QGRAM_LENGTH_POOLS_QUALITY: {
                auto tmp = new LengthPoolsAmpliconStorage<ArrayQgramQualityAmpliconCollection>();

                std::vector<lenSeqs_t> lengths = ds.get_all_lengths();
                auto pooling_threshold = static_cast<lenSeqs_t>(std::max(config.main_threshold, config.refinement_threshold));

                lenSeqs_t last_len = ds.get_min_length();
                numSeqs_t pool_id = 0;

                numSeqs_t cnt = 0;
                unsigned long long total_len_sequences = 0;
                unsigned long long total_len_headers = 0;

                for (auto len : lengths) {

                    if ((last_len + pooling_threshold) < len) { // new pool

                        pool_id++;
                        tmp->pools_.emplace_back(cnt, total_len_headers, total_len_sequences);

                        cnt = 0;
                        total_len_sequences = 0;
                        total_len_headers = 0;

                    }

                    cnt += ds.get_num_per_length(len);
                    total_len_sequences += (len + 1) * ds.get_num_per_length(len);
                    total_len_headers += ds.get_total_length_headers_per_length(len) + ds.get_num_per_length(len);

                    tmp->pool_map_[len] = pool_id;
                    last_len = len;

                }

                tmp->pools_.emplace_back(cnt, total_len_headers, total_len_sequences);
                res = tmp;
                break;
            }

            case AS_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown AmpliconStorage option." << std::endl;

        }

        return res;

    }

    /*
     * The proper initialisation of the more involved data structures considered in the investigation
     * of space efficiency require additional pool-wise information.
     */
    struct PoolsInformation {

        numSeqs_t num_pools = 0; // number of pools
        std::vector<numSeqs_t> size; // size of each pool
        std::vector<unsigned long long> total_len_headers; // total length of headers in each pool
        std::vector<unsigned long long> total_len_sequences; // total length of sequences in each pool
        std::map<lenSeqs_t, numSeqs_t> pool_map; // mapping between sequence length and pool
        std::vector<std::array<size_t, 256>> char_counts; // distribution of characters in headers per pool
        std::vector<std::map<ularge_t, ularge_t>> abundance_counts; // distribution of abundance values per pool

        /*
         * Add the information of next pool.
         */
        void add_pool(const numSeqs_t s, const numSeqs_t tlh, const numSeqs_t tls,
                      const std::array<size_t, 256>& cc, const std::map<ularge_t, ularge_t>& ac) {

            num_pools++;
            size.push_back(s);
            total_len_headers.push_back(tlh);
            total_len_sequences.push_back(tls);
            char_counts.push_back(cc);
            abundance_counts.push_back(ac);

        }

        /*
         * Extend the length-pool mapping.
         */
        void add_length(lenSeqs_t len, numSeqs_t pi) {
            pool_map[len] = pi;
        }

        /*
         * Determine all lengths that belong to a pool.
         */
        std::vector<lenSeqs_t> lengths_in_pool(numSeqs_t pi) {

            std::vector<lenSeqs_t> lengths;
            for (auto& kv : pool_map) if (kv.second == pi) lengths.push_back(kv.first);
            return lengths;

        }

    };

    /*
     * Determine the information on pools from the data statistics and the configuration.
     */
    PoolsInformation pool(const DataStatistics<>& ds, const Configuration& config) {

        PoolsInformation pools;
        std::vector<lenSeqs_t> lengths = ds.get_all_lengths();
        auto pooling_threshold = static_cast<lenSeqs_t>(std::max(config.main_threshold, config.refinement_threshold));

        lenSeqs_t last_len = ds.get_min_length();
        numSeqs_t pool_id = 0;

        numSeqs_t cnt = 0;
        unsigned long long total_len_sequences = 0;
        unsigned long long total_len_headers = 0;
        std::array<size_t, 256> char_counts = {};
        std::map<ularge_t, ularge_t> abundance_counts;

        for (auto len : lengths) {

            if ((last_len + pooling_threshold) < len) { // new pool

                pool_id++;
                pools.add_pool(cnt, total_len_headers, total_len_sequences, char_counts, abundance_counts);

                cnt = 0;
                total_len_sequences = 0;
                total_len_headers = 0;
                char_counts = {};
                abundance_counts.clear();

            }

            cnt += ds.get_num_per_length(len);
            total_len_sequences += (len + 1) * ds.get_num_per_length(len);
            total_len_headers += ds.get_total_length_headers_per_length(len) + ds.get_num_per_length(len);
            auto& counts = ds.get_counts_per_length(len);
            for (auto i = 0; i < 256; i++) char_counts[i] += counts[i];
            auto& ab_counts = ds.get_abundances_per_length(len);
            for (auto& kv : ab_counts) abundance_counts[kv.first] += kv.second;

            pools.add_length(len, pool_id);
            last_len = len;

        }

        pools.add_pool(cnt, total_len_headers, total_len_sequences, char_counts, abundance_counts);

        return pools;

    }

    AmpliconStorage* AmpliconStorageFactory::create(const AmpliconStorageOption opt, const AmpliconCollectionOption ac_opt,
                                                    const DataStatistics<>& ds, const SpaceLevenshteinConfiguration& config) {

        AmpliconStorage* res = nullptr;

        switch (opt) {

            // Reference:
            // Create multiple pools when there are sufficient gaps in the length distribution,
            // across which amplicons cannot be similar.
            // Each pool is initialised with a sufficient capacity for the amplicons (including their q-gram profiles)
            // expected to be inserted.
            case AS_PREPARED_QGRAM_LENGTH_POOLS: {
                auto tmp = new LengthPoolsAmpliconStorage<ArrayQgramAmpliconCollection>();

                std::vector<lenSeqs_t> lengths = ds.get_all_lengths();
                auto pooling_threshold = static_cast<lenSeqs_t>(std::max(config.main_threshold, config.refinement_threshold));

                lenSeqs_t last_len = ds.get_min_length();
                numSeqs_t pool_id = 0;

                numSeqs_t cnt = 0;
                unsigned long long total_len_sequences = 0;
                unsigned long long total_len_headers = 0;

                for (auto len : lengths) {

                    if ((last_len + pooling_threshold) < len) { // new pool

                        pool_id++;
                        tmp->pools_.emplace_back(cnt, total_len_headers, total_len_sequences);

                        cnt = 0;
                        total_len_sequences = 0;
                        total_len_headers = 0;

                    }

                    cnt += ds.get_num_per_length(len);
                    total_len_sequences += (len + 1) * ds.get_num_per_length(len);
                    total_len_headers += ds.get_total_length_headers_per_length(len) + ds.get_num_per_length(len);

                    tmp->pool_map_[len] = pool_id;
                    last_len = len;

                }

                tmp->pools_.emplace_back(cnt, total_len_headers, total_len_sequences);
                res = tmp;
                break;
            }

            // Space-effient variants of amplicon storage:
            // Create multiple pools as for AS_PREPARED_QGRAM_LENGTH_POOLS.
            // Build components for headers, sequences, lengths, abundances and q-gram profiles based on the configuration.
            case AS_SPACE_FLEXIBLE: {
                PoolsInformation pools = pool(ds, config);

                auto tmp = new LengthPoolsAmpliconStorage<FlexibleAmpliconCollection>();
                for (auto& kv : pools.pool_map) tmp->pool_map_[kv.first] = kv.second;

                switch (ac_opt) {

                    case AC_SPACE_REFERENCE: {
                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            ac.set_header_component(new CollectiveStringArray(pools.size[i], pools.total_len_headers[i]));
                            ac.set_sequence_component(new CollectiveStringArray(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new SimpleArray<numSeqs_t>(pools.size[i]));
                            ac.set_length_component(new SimpleArray<lenSeqs_t>(pools.size[i]));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_ID_PREF: {
                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            ac.set_header_component(new PrefixStringArray(config.prefix_len, pools.size[i], pools.total_len_headers[i]));
                            ac.set_sequence_component(new CollectiveStringArray(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new SimpleArray<numSeqs_t>(pools.size[i]));
                            ac.set_length_component(new SimpleArray<lenSeqs_t>(pools.size[i]));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_ID_PREF_DACS: {
                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            ac.set_header_component(new PrefixDacsStringArray(config.prefix_len, pools.size[i], pools.total_len_headers[i], config.chunk_lengths));
                            ac.set_sequence_component(new CollectiveStringArray(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new SimpleArray<numSeqs_t>(pools.size[i]));
                            ac.set_length_component(new SimpleArray<lenSeqs_t>(pools.size[i]));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_ID_HUFF: {
                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            ac.set_header_component(new HuffmanStringArray<size_t>(pools.size[i], pools.char_counts[i]));
                            ac.set_sequence_component(new CollectiveStringArray(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new SimpleArray<numSeqs_t>(pools.size[i]));
                            ac.set_length_component(new SimpleArray<lenSeqs_t>(pools.size[i]));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_ID_HUFF_DACS: {
                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            ac.set_header_component(new HuffmanDacsStringArray<size_t>(pools.size[i], pools.char_counts[i], config.chunk_lengths));
                            ac.set_sequence_component(new CollectiveStringArray(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new SimpleArray<numSeqs_t>(pools.size[i]));
                            ac.set_length_component(new SimpleArray<lenSeqs_t>(pools.size[i]));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_SEQ: {
                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            ac.set_header_component(new CollectiveStringArray(pools.size[i], pools.total_len_headers[i]));
                            ac.set_sequence_component(new TwoBitStringArray<size_t>(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new SimpleArray<numSeqs_t>(pools.size[i]));
                            ac.set_length_component(new SimpleArray<lenSeqs_t>(pools.size[i]));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_SEQ_DACS: {
                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            ac.set_header_component(new CollectiveStringArray(pools.size[i], pools.total_len_headers[i]));
                            ac.set_sequence_component(new TwoBitDacsStringArray<size_t>(pools.size[i], pools.total_len_sequences[i], config.chunk_lengths));
                            ac.set_abundance_component(new SimpleArray<numSeqs_t>(pools.size[i]));
                            ac.set_length_component(new SimpleArray<lenSeqs_t>(pools.size[i]));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_AB: {
                        typedef Dacs<numSeqs_t, size_t> dacs_t;
                        typedef DacsBuilder<numSeqs_t, size_t> builder_t;

                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            ac.set_header_component(new CollectiveStringArray(pools.size[i], pools.total_len_headers[i]));
                            ac.set_sequence_component(new CollectiveStringArray(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new DacsArray<dacs_t, builder_t, numSeqs_t>(pools.abundance_counts[i]));
                            ac.set_length_component(new SimpleArray<lenSeqs_t>(pools.size[i]));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_AB_BITS: {
                        typedef Dacs<numSeqs_t, size_t> dacs_t;
                        typedef DacsBuilder<numSeqs_t, size_t> builder_t;

                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            ac.set_header_component(new CollectiveStringArray(pools.size[i], pools.total_len_headers[i]));
                            ac.set_sequence_component(new CollectiveStringArray(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new DacsBitsArray<dacs_t, builder_t, numSeqs_t, size_t>(pools.abundance_counts[i]));
                            ac.set_length_component(new SimpleArray<lenSeqs_t>(pools.size[i]));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_LEN: {
                        typedef Dacs<numSeqs_t, size_t> dacs_t;
                        typedef DacsBuilder<numSeqs_t, size_t> builder_t;

                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            std::map<ularge_t, ularge_t> counts;
                            for (auto len : pools.lengths_in_pool(i)) counts[len] = ds.get_num_per_length(len);

                            ac.set_header_component(new CollectiveStringArray(pools.size[i], pools.total_len_headers[i]));
                            ac.set_sequence_component(new CollectiveStringArray(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new SimpleArray<numSeqs_t>(pools.size[i]));
                            ac.set_length_component(new DacsArray<dacs_t, builder_t, lenSeqs_t>(counts));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_LEN_OFFSET: {
                        typedef Dacs<numSeqs_t, size_t> dacs_t;
                        typedef DacsBuilder<numSeqs_t, size_t> builder_t;

                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            std::map<ularge_t, ularge_t> counts;
                            for (auto len : pools.lengths_in_pool(i)) counts[len] = ds.get_num_per_length(len);

                            ac.set_header_component(new CollectiveStringArray(pools.size[i], pools.total_len_headers[i]));
                            ac.set_sequence_component(new CollectiveStringArray(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new SimpleArray<numSeqs_t>(pools.size[i]));
                            ac.set_length_component(new DacsOffsetArray<dacs_t, builder_t, lenSeqs_t>(counts));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_AB_SPLIT: {
                        typedef SplitDacs<numSeqs_t, size_t> dacs_t;
                        typedef SplitDacsDistributionBuilder<numSeqs_t, size_t> builder_t;

                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            ac.set_header_component(new CollectiveStringArray(pools.size[i], pools.total_len_headers[i]));
                            ac.set_sequence_component(new CollectiveStringArray(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new DacsArray<dacs_t, builder_t, numSeqs_t>(pools.abundance_counts[i]));
                            ac.set_length_component(new SimpleArray<lenSeqs_t>(pools.size[i]));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_AB_BITS_SPLIT: {
                        typedef SplitDacs<numSeqs_t, size_t> dacs_t;
                        typedef SplitDacsDistributionBuilder<numSeqs_t, size_t> builder_t;

                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            ac.set_header_component(new CollectiveStringArray(pools.size[i], pools.total_len_headers[i]));
                            ac.set_sequence_component(new CollectiveStringArray(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new DacsBitsArray<dacs_t, builder_t, numSeqs_t, size_t>(pools.abundance_counts[i]));
                            ac.set_length_component(new SimpleArray<lenSeqs_t>(pools.size[i]));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_LEN_SPLIT: {
                        typedef SplitDacs<numSeqs_t, size_t> dacs_t;
                        typedef SplitDacsDistributionBuilder<numSeqs_t, size_t> builder_t;

                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            std::map<ularge_t, ularge_t> counts;
                            for (auto len : pools.lengths_in_pool(i)) counts[len] = ds.get_num_per_length(len);

                            ac.set_header_component(new CollectiveStringArray(pools.size[i], pools.total_len_headers[i]));
                            ac.set_sequence_component(new CollectiveStringArray(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new SimpleArray<numSeqs_t>(pools.size[i]));
                            ac.set_length_component(new DacsArray<dacs_t, builder_t, lenSeqs_t>(counts));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_LEN_OFFSET_SPLIT: {
                        typedef SplitDacs<numSeqs_t, size_t> dacs_t;
                        typedef SplitDacsDistributionBuilder<numSeqs_t, size_t> builder_t;

                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            std::map<ularge_t, ularge_t> counts;
                            for (auto len : pools.lengths_in_pool(i)) counts[len] = ds.get_num_per_length(len);

                            ac.set_header_component(new CollectiveStringArray(pools.size[i], pools.total_len_headers[i]));
                            ac.set_sequence_component(new CollectiveStringArray(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new SimpleArray<numSeqs_t>(pools.size[i]));
                            ac.set_length_component(new DacsOffsetArray<dacs_t, builder_t, lenSeqs_t>(counts));
                            ac.set_qgram_component(new QgramsArray(pools.size[i]));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_QGRAM_DACS: {
                        typedef SplitDacs<numSeqs_t, size_t> dacs_t;
                        typedef SplitDacsDistributionBuilder<numSeqs_t, size_t> builder_t;

                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            std::map<ularge_t, ularge_t> counts;
                            for (auto len : pools.lengths_in_pool(i)) counts[len] = ds.get_num_per_length(len);

                            ac.set_header_component(new CollectiveStringArray(pools.size[i], pools.total_len_headers[i]));
                            ac.set_sequence_component(new CollectiveStringArray(pools.size[i], pools.total_len_sequences[i]));
                            ac.set_abundance_component(new SimpleArray<numSeqs_t>(pools.size[i]));
                            ac.set_length_component(new SimpleArray<lenSeqs_t>(pools.size[i]));
                            ac.set_qgram_component(new QgramsDacsArray(pools.size[i], config.chunk_lengths));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_SPACE_FULL: {
                        typedef Dacs<numSeqs_t, size_t> dacs_t;
                        typedef DacsBuilder<numSeqs_t, size_t> builder_t;

                        for (auto i = 0; i < pools.num_pools; i++) {

                            tmp->pools_.emplace_back(pools.size[i]);
                            auto& ac = tmp->pools_.back();

                            std::map<ularge_t, ularge_t> counts;
                            for (auto len : pools.lengths_in_pool(i)) counts[len] = ds.get_num_per_length(len);

                            ac.set_header_component(new PrefixDacsStringArray(config.prefix_len, pools.size[i], pools.total_len_headers[i], config.chunk_lengths));
                            ac.set_sequence_component(new TwoBitDacsStringArray<size_t>(pools.size[i], pools.total_len_sequences[i], config.chunk_lengths));
                            ac.set_abundance_component(new DacsBitsArray<dacs_t, builder_t, numSeqs_t, size_t>(pools.abundance_counts[i]));
                            ac.set_length_component(new DacsOffsetArray<dacs_t, builder_t, lenSeqs_t>(counts));
                            ac.set_qgram_component(new QgramsDacsArray(pools.size[i], config.chunk_lengths));

                        }
                        res = tmp;
                        break;
                    }

                    case AC_UNKNOWN: // fallthrough to default for error handling
                    default:
                        std::cerr << "ERROR: Unknown AmpliconCollection option." << std::endl;

                }
                break;
            }

            case AS_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown AmpliconStorage option." << std::endl;

        }

        return res;

    }

    FeatureAmpliconStorage* AmpliconStorageFactory::create(const AmpliconStorageOption opt, const FeatureBuilderOption fb_opt,
            const AmpliconCollectionOption ac_opt, const DataStatistics<>& ds, const Configuration& config) {

        FeatureBuilder* fb = nullptr;
        auto& misc = config.misc;

        switch (fb_opt) {

            case FB_WCV: {

                short k = static_cast<short>(std::stoi(config.misc.find("wcv_word_length")->second)); // existence ensured by AlignmentFreeConfiguration
                fb = new WordCompBuilder(config.alphabet, k);

                break;
            }

            case FB_CPF: {

                fb = new CpfBuilder(ds.get_max_length());

                break;
            }

            case FB_DET: {

                // existence of either parameter ensured by AlignmentFreeConfiguration
                if (config.misc.find("det_alpha") != config.misc.end()) {
                    feat_t alpha = std::stof(config.misc.find("det_alpha")->second);
                    fb = new DetBuilder(config.alphabet, alpha, ds.get_max_length());
                }
                if (config.misc.find("det_pref_dist") != config.misc.end()) {
                    lenSeqs_t pref_dist = std::stoul(config.misc.find("det_pref_dist")->second);
                    fb = new DetBuilder(config.alphabet, pref_dist, ds.get_max_length());
                }

                break;
            }

            case FB_BBC: {

                short max_dist = static_cast<short>(std::stoi(config.misc.find("bbc_max_dist")->second)); // existence ensured by AlignmentFreeConfiguration
                fb = new BbcBuilder(config.alphabet, max_dist);

                break;
            }

            case FB_2DW: {

                feat_t alpha = std::stof(config.misc.find("2dw_alpha")->second); // existence ensured by AlignmentFreeConfiguration
                feat_t beta = std::stof(config.misc.find("2dw_beta")->second); // existence ensured by AlignmentFreeConfiguration
                fb = new TwoDimWalkBuilder(alpha, beta);

                break;
            }

            case FB_CGR: {

                fb = new CgrBuilder();

                break;
            }

            case FB_MCGR: {

                fb = new MultiCgrBuilder();

                break;
            }

            case FB_3DCGR: {

                fb = new ThreeDimCgrBuilder();

                break;
            }

            case FB_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown FeatureBuilder option." << std::endl;

        }

        FeatureAmpliconStorage* res = nullptr;

        if (fb != nullptr) {

            switch (opt) {

                // Create a single pool for all amplicons.
                case AS_SIMPLE_FEATURES: {

                    switch (ac_opt) {

                        case AC_ARRAY: {
                            std::vector<lenSeqs_t> lengths = ds.get_all_lengths();

                            numSeqs_t cnt = 0;
                            unsigned long long total_len_sequences = 0;
                            unsigned long long total_len_headers = 0;

                            for (auto len : lengths) {

                                cnt += ds.get_num_per_length(len);
                                total_len_sequences += (len + 1) * ds.get_num_per_length(len);
                                total_len_headers += ds.get_total_length_headers_per_length(len) + ds.get_num_per_length(len);

                            }

                            ArrayFeatVecAmpliconCollection pool(fb, cnt, total_len_headers, total_len_sequences);
                            res = new SimpleFeatureAmpliconStorage<ArrayFeatVecAmpliconCollection>(pool);
                            break;
                        }

                        case AC_SIMPLE: // fallthrough
                        default:
                            SimpleFeatVecAmpliconCollection pool(fb);
                            res = new SimpleFeatureAmpliconStorage<SimpleFeatVecAmpliconCollection>(pool);

                    }
                    break;

                }

                case AS_UNKNOWN: // fallthrough to default for error handling
                default:
                    std::cerr << "ERROR: Unknown AmpliconStorage option." << std::endl;

            }

        }

        return res;

    }


    Distance* DistanceFactory::create(const DistanceOption opt, const AmpliconStorage& amplicon_storage,
            const dist_t threshold, const Configuration& config) {

        Distance* res = nullptr;

        switch (opt) {

            case DT_BOUNDED_LEVENSHTEIN: {
                res = new BoundedLevenshteinDistance(amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold));
                break;
            }

            case DT_BOUNDED_SCORE_LEVENSHTEIN: {
                res = new BoundedScoreLevenshteinDistance(amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold),
                        config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
                break;
            }

            case DT_BOUNDED_SCORE: {
                res = new BoundedScoreDistance(amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold),
                        config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
                break;
            }

            case DT_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown Distance option." << std::endl;

        }

        return res;

    }


    Distance* DistanceFactory::create(const DistanceOption opt, const AmpliconStorage& amplicon_storage,
            const dist_t threshold, const LevenshteinConfiguration& config) {

        Distance* res = nullptr;

        switch (opt) {

            case DT_BOUNDED_LEVENSHTEIN: {
                res = new BoundedLevenshteinDistance(amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold));
                break;
            }

            case DT_BOUNDED_SCORE_LEVENSHTEIN: {
                res = new BoundedScoreLevenshteinDistance(amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold),
                        config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
                break;
            }

            case DT_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown Distance option for Levenshtein (lev) mode." << std::endl;

        }

        if (res != nullptr && config.use_qgrams) { // wraps the distance computation by a check of the q-gram distance
            res = new QgramBoundedLevenshteinDistance(res, static_cast<lenSeqs_t>(threshold));
        }

        return res;

    }

    Distance* DistanceFactory::create(const DistanceOption opt, const AmpliconStorage& amplicon_storage,
            const dist_t threshold, const AlignmentScoreConfiguration& config) {

        Distance* res = nullptr;

        switch (opt) {

            case DT_BOUNDED_SCORE: {
                res = new BoundedScoreDistance(amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold),
                        config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
                break;
            }

            case DT_BOUNDED_BANDED_SCORE: {
                long long bands_per_side = config.bands_per_side;
                if (bands_per_side == -1) {

                    ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
                    bands_per_side = static_cast<lenSeqs_t>(threshold) / std::min({sf.penalty_mismatch, sf.penalty_open, sf.penalty_extend});

                }

                res = new BoundedBandedScoreDistance(amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold),
                        bands_per_side, config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);

                break;
            }

            case DT_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown Distance option for alignment-score (as) mode." << std::endl;

        }

        if (res != nullptr && config.use_qgrams) { // wraps the distance computation by a check of the q-gram distance
            res = new QgramBoundedScoreDistance(res, static_cast<lenSeqs_t>(threshold), config.match_reward,
                    config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
        }

        return res;

    }


    Distance* DistanceFactory::create(const DistanceOption opt, const AmpliconStorage& amplicon_storage,
            const dist_t threshold, const QualityAlignmentScoreConfiguration& config) {

        Distance* res = nullptr;

        long long bands_per_side = config.bands_per_side;
        if (bands_per_side == -1) {

            ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
            bands_per_side = static_cast<lenSeqs_t>(threshold) / std::min({sf.penalty_mismatch, sf.penalty_open, sf.penalty_extend});

        }

        switch (opt) {

            case DT_QWS_CLEMENT: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityAlignmentScore, ClementScores>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityAlignmentScore, ClementScoresWithMatches>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QWS_MALDE_A: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityAlignmentScore, MaldeScoresA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityAlignmentScore, MaldeScoresWithMatchesA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QWS_MALDE_B: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityAlignmentScore, MaldeScoresB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityAlignmentScore, MaldeScoresWithMatchesB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QWS_MALDE_C: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityAlignmentScore, MaldeScoresC>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityAlignmentScore, MaldeScoresWithMatchesC>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QWS_FRITH: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityAlignmentScore, FrithScores>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityAlignmentScore, FrithScoresWithMatches>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QWS_KIM_A: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityAlignmentScore, KimScoresA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityAlignmentScore, KimScoresWithMatchesA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QWS_KIM_B: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityAlignmentScore, KimScoresB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityAlignmentScore, KimScoresWithMatchesB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QWS_CONVERGE_A: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityAlignmentScore, ConvergeScoresA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityAlignmentScore, ConvergeScoresWithMatchesA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QWS_CONVERGE_B: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityAlignmentScore, ConvergeScoresB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityAlignmentScore, ConvergeScoresWithMatchesB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }


            case DT_QWBS_CLEMENT: {
                res = (config.unweighted_matches) ?
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, ClementScores>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config) :
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, ClementScoresWithMatches>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config);
                break;
            }

            case DT_QWBS_MALDE_A: {
                res = (config.unweighted_matches) ?
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, MaldeScoresA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config) :
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, MaldeScoresWithMatchesA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config);
                break;
            }

            case DT_QWBS_MALDE_B: {
                res = (config.unweighted_matches) ?
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, MaldeScoresB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config) :
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, MaldeScoresWithMatchesB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config);
                break;
            }

            case DT_QWBS_MALDE_C: {
                res = (config.unweighted_matches) ?
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, MaldeScoresC>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config) :
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, MaldeScoresWithMatchesC>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config);
                break;
            }

            case DT_QWBS_FRITH: {
                res = (config.unweighted_matches) ?
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, FrithScores>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config) :
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, FrithScoresWithMatches>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config);
                break;
            }

            case DT_QWBS_KIM_A: {
                res = (config.unweighted_matches) ?
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, KimScoresA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config) :
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, KimScoresWithMatchesA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config);
                break;
            }

            case DT_QWBS_KIM_B: {
                res = (config.unweighted_matches) ?
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, KimScoresB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config) :
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, KimScoresWithMatchesB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config);
                break;
            }

            case DT_QWBS_CONVERGE_A: {
                res = (config.unweighted_matches) ?
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, ConvergeScoresA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config) :
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, ConvergeScoresWithMatchesA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config);
                break;
            }

            case DT_QWBS_CONVERGE_B: {
                res = (config.unweighted_matches) ?
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, ConvergeScoresB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config) :
                      get_bqw_distance<BoundedBandedQualityAlignmentScore, ConvergeScoresWithMatchesB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), bands_per_side, config);
                break;
            }

            case DT_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown Distance option for quality alignment-score (qas) mode." << std::endl;

        }

        if (res != nullptr && config.use_qgrams) { // wraps the distance computation by a check of the q-gram distance
            res = new QgramBoundedScoreDistance(res, static_cast<lenSeqs_t>(threshold), config.match_reward,
                                                config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
        }

        return res;

    }

    Distance* DistanceFactory::create(const DistanceOption opt, const AmpliconStorage& amplicon_storage, const dist_t threshold,
            const QualityLevenshteinConfiguration& config) {

        Distance* res = nullptr;

        switch (opt) {

            case DT_QW_CLEMENT: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityLevenshteinDistance, ClementScores>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityLevenshteinDistance, ClementScoresWithMatches>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QW_MALDE_A: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityLevenshteinDistance, MaldeScoresA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityLevenshteinDistance, MaldeScoresWithMatchesA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QW_MALDE_B: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityLevenshteinDistance, MaldeScoresB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityLevenshteinDistance, MaldeScoresWithMatchesB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QW_MALDE_C: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityLevenshteinDistance, MaldeScoresC>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityLevenshteinDistance, MaldeScoresWithMatchesC>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QW_FRITH: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityLevenshteinDistance, FrithScores>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityLevenshteinDistance, FrithScoresWithMatches>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QW_KIM_A: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityLevenshteinDistance, KimScoresA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityLevenshteinDistance, KimScoresWithMatchesA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QW_KIM_B: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityLevenshteinDistance, KimScoresB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityLevenshteinDistance, KimScoresWithMatchesB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QW_CONVERGE_A: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityLevenshteinDistance, ConvergeScoresA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityLevenshteinDistance, ConvergeScoresWithMatchesA>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_QW_CONVERGE_B: {
                res = (config.unweighted_matches) ?
                      get_qw_distance<BoundedQualityLevenshteinDistance, ConvergeScoresB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config) :
                      get_qw_distance<BoundedQualityLevenshteinDistance, ConvergeScoresWithMatchesB>(
                              amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold), config);
                break;
            }

            case DT_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown Distance option for quality Levenshtein (qlev) mode." << std::endl;

        }

        if (res != nullptr && config.use_qgrams) { // wraps the distance computation by a check of the q-gram distance
            res = new QgramBoundedLevenshteinDistance(res, static_cast<lenSeqs_t>(threshold));
        }

        return res;

    }

    FeatureDistance* DistanceFactory::create(const DistanceOption opt, const AmpliconStorage& amplicon_storage,
                                             const dist_t threshold, const AlignmentFreeConfiguration& config) {

        FeatureDistance* res = nullptr;

        switch (opt) {

            case DT_EUCLIDEAN: {
                res = new EuclideanDistance();
                break;
            }

            case DT_MANHATTAN: {
                res = new ManhattanDistance();
                break;
            }

            case DT_COSINE: {
                res = new CosineDistance();
                break;
            }

            case DT_PEARSON: {
                res = new PearsonDistance();
                break;
            }

            case DT_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown Distance option." << std::endl;

        }

        return res;

    }


    Distance* DistanceFactory::create(const DistanceOption opt, const AmpliconStorage& amplicon_storage,
                                      const dist_t threshold, const SpaceLevenshteinConfiguration& config) {

        Distance* res = nullptr;

        switch (opt) {

            case DT_BOUNDED_LEVENSHTEIN: {
                res = new BoundedLevenshteinDistance(amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold));
                break;
            }

            case DT_BOUNDED_SCORE_LEVENSHTEIN: {
                res = new BoundedScoreLevenshteinDistance(amplicon_storage.max_length(), static_cast<lenSeqs_t>(threshold),
                                                          config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
                break;
            }

            case DT_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown Distance option for space-efficiency (space) mode." << std::endl;

        }

        if (res != nullptr && config.use_qgrams) { // wraps the distance computation by a check of the q-gram distance
            res = new QgramBoundedLevenshteinDistance(res, static_cast<lenSeqs_t>(threshold));
        }

        return res;

    }


    Preprocessor* PreprocessorFactory::create(const PreprocessorOption opt, const QualityEncoding<>& qe,
            const Configuration& config) {

        Preprocessor* res = nullptr;

        switch (opt) {
            case PP_FASTA:
                res = new FastaPreprocessor();
                break;

            case PP_FASTQ: {
                auto iter = config.misc.find("score_derep");
                std::string derep_opt = (iter != config.misc.end()) ? iter->second : "mean";

                if (derep_opt == "mean") {
                    res = new FastqPreprocessor<MeanReadRepresentation>(qe);
                } else if (derep_opt == "max") {
                    res = new FastqPreprocessor<MaxReadRepresentation>(qe);
                } else if (derep_opt == "mode") {
                    res = new FastqPreprocessor<ModeReadRepresentation>(qe);
                } else if (derep_opt == "median") {
                    res = new FastqPreprocessor<MedianReadRepresentation>(qe);
                } else if (derep_opt == "product") {
                    res = new FastqPreprocessor<ProductReadRepresentation>(qe);
                } else {
                    std::cerr << "ERROR: Unknown 'score_derep' option. Resorting to 'mean'." << std::endl;
                    res = new FastqPreprocessor<MeanReadRepresentation>(qe);
                }
                break;
            }

            case PP_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown preprocessor type." << std::endl;

        }

        return res;

    }

    Preprocessor* PreprocessorFactory::create(const PreprocessorOption opt_prep, const QualityEncodingOption opt_qe,
            const Configuration& config) {

        Preprocessor* res = nullptr;

        switch (opt_prep) {
            case PP_FASTA:
                res = new FastaPreprocessor();
                break;

            case PP_FASTQ: {
                auto iter = config.misc.find("score_derep");
                std::string derep_opt = (iter != config.misc.end()) ? iter->second : "mean";

                if (derep_opt == "mean") {
                    res = new FastqPreprocessor<MeanReadRepresentation>(opt_qe);
                } else if (derep_opt == "max") {
                    res = new FastqPreprocessor<MaxReadRepresentation>(opt_qe);
                } else if (derep_opt == "mode") {
                    res = new FastqPreprocessor<ModeReadRepresentation>(opt_qe);
                } else if (derep_opt == "median") {
                    res = new FastqPreprocessor<MedianReadRepresentation>(opt_qe);
                } else if (derep_opt == "product") {
                    res = new FastqPreprocessor<ProductReadRepresentation>(opt_qe);
                } else {
                    std::cerr << "ERROR: Unknown 'score_derep' option. Resorting to 'mean'." << std::endl;
                    res = new FastqPreprocessor<MeanReadRepresentation>(opt_qe);
                }
                break;
            }

            case PP_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown preprocessor type." << std::endl;

        }

        return res;

    }



    Clusterer* ClustererFactory::create(const ClustererOption opt, const Configuration& config) {

        Clusterer* res = nullptr;

        switch (opt) {

            case CL_IDLE: {
                res = new IdleClusterer();
                break;
            }

            case CL_CLASSIC_SWARMER: {
                res = new ClassicSwarmer();
                break;
            }

            case CL_DEREPLICATOR: {
                res = new Dereplicator();
                break;
            }

            case CL_DADA2: {
                res = new Dada2Clusterer(config);
                break;
            }

            case CL_CONS_CLASSIC: {
                res = new ConsistentClassicSwarmer(config);
                break;
            }

            case CL_CONS_SWARMER: {
                res = new ConsistencySwarmer(config);
                break;
            }

            case CL_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown Clusterer option." << std::endl;

        }

        return res;

    }



    ClusterRefiner* ClusterRefinerFactory::create(const ClusterRefinerOption opt, const Configuration& config) {

        ClusterRefiner* res = nullptr;

        switch (opt) {

            case CR_IDLE: {
                res = new IdleRefiner();
                break;
            }

            case CR_CLASSIC_FASTIDIOUS: {
                res = new FastidiousRefiner();
                break;
            }

            case CR_ITERATIVE_FASTIDIOUS: {
                res = new IterativeFastidiousRefiner();
                break;
            }

            case CR_CONS_LSA: {
                res = new LightSwarmAppender(config);
                break;
            }

            case CR_CONS_LSR: {
                res = new LightSwarmResolver(config);
                break;
            }

            case CR_CONS_LSS: {
                res = new LightSwarmShuffler(config);
                break;
            }

            case CR_CONS_SS: {
                res = new SwarmShuffler(config);
                break;
            }

            case CR_UNKNOWN:
            default:
                std::cerr << "ERROR: Unknown ClusterRefiner option." << std::endl;

        }

        return res;

    }



    OutputGenerator* OutputGeneratorFactory::create(const OutputGeneratorOption opt, const Configuration& config) {

        OutputGenerator* res = nullptr;

        switch (opt) {

            case OG_IDLE: {
                res = new IdleOutputGenerator();
                break;
            }

            case OG_CLASSIC: {
                res = new ClassicOutputGenerator();
                break;
            }

            case OG_DEREPLICATION: {
                res = new DereplicationOutputGenerator();
                break;
            }

            case OG_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown OutputGenerator option." << std::endl;

        }

        return res;

    }



    SwarmStorage* SwarmStorageFactory::create(const SwarmStorageOption opt, const AmpliconStorage& amplicon_storage,
            const Configuration& config) {

        SwarmStorage* res = nullptr;

        switch (opt) {

            case SS_SIMPLE_PER_POOL: {
                res = new SimpleSwarmStorage<SimpleSwarms<SimpleSwarm>>(amplicon_storage.num_pools());
                break;
            }

            case SS_SIMPLE_COMB_PER_POOL: {
                res = new SimpleSwarmStorage<SimpleSwarms<SimpleCombinedSwarm>>(amplicon_storage.num_pools());
                break;
            }

            case SS_FORWARDING_PER_POOL: {
                res = new SimpleSwarmStorage<SharingSplitSwarms<ForwardingSwarm>>(amplicon_storage.num_pools());
                break;
            }

            case SS_FORWARDING_COMB_PER_POOL: {
                res = new SimpleSwarmStorage<SharingCombinedSwarms<ForwardingSwarm>>(amplicon_storage.num_pools());
                break;
            }

            case SS_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown SwarmStorage option." << std::endl;

        }

        return res;

    }



    SwarmStorage* SwarmStorageFactory::create(const SwarmStorageOption opt, const AmpliconStorage& amplicon_storage,
            const SpaceLevenshteinConfiguration& config) {

        SwarmStorage* res = nullptr;

        switch (opt) {

            case SS_FORWARDING_DACS_PER_POOL: {
                res = new DacsSwarmStorage<SharingSplitDacsSwarms>(amplicon_storage, config);
                break;
            }

            case SS_FORWARDING_FULL_DACS_PER_POOL: {
                res = new DacsSwarmStorage<SharingSplitFullDacsSwarms>(amplicon_storage, config);
                break;
            }

            default:
                res = create(opt, amplicon_storage, static_cast<const Configuration&>(config));

        }

        return res;

    }



    AuxiliaryData* AuxiliaryDataFactory::create(const AuxiliaryDataOption opt, const AmpliconStorage& amplicon_storage,
            const numSeqs_t pool_id, const dist_t threshold, const Configuration& config) {

        AuxiliaryData* res = nullptr;

        switch (opt) {

            // Auxiliary data structures providing the bare minimum to execute the clustering phase.
            // No filtering of potential partners.
            case AD_NAIVE_AUXILIARY: {
                res = new SimpleAuxiliaryData(amplicon_storage, pool_id, threshold, config);
                break;
            }

            case AD_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown AuxiliaryData option." << std::endl;

        }

        return res;

    }

    AuxiliaryData* AuxiliaryDataFactory::create(const AuxiliaryDataOption opt, const AmpliconStorage& amplicon_storage,
            const numSeqs_t pool_id, const dist_t threshold, const lenSeqs_t num_extra_segments, const Configuration& config) {

        AuxiliaryData* res = nullptr;

        switch (opt) {

            // Auxiliary data structures employing a (one-way) segment filter to efficiently search
            // for (not yet swarmed) partners of an amplicon.
            case AD_SEGMENT_FILTER: {
                res = new SegmentFilterAuxiliaryData(amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                        num_extra_segments, config);
                break;
            }

            // Auxiliary data structures employing a (two-way) segment filter to efficiently search
            // for (not yet swarmed) partners of an amplicon.
            case AD_2W_SEGMENT_FILTER: {
                res = new TwoWaySegmentFilterAuxiliaryData(amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                        num_extra_segments, config);
                break;
            }

            // Auxiliary data structures employing a (one-way) segment filter to efficiently search
            // for (not yet swarmed) partners of an amplicon.
            // The set-up process of the segment filter is adapted in order to determine the number of segments
            // needed for the filter based on the distance threshold and the edit-operation costs.
            case AD_SCORE_SEGMENT_FILTER: {
                ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
                lenSeqs_t num_threshold_segments = static_cast<lenSeqs_t>(threshold) / std::min({sf.penalty_mismatch, sf.penalty_open, sf.penalty_extend});

                res = new ScoreSegmentFilterAuxiliaryData(amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                        num_threshold_segments, num_extra_segments, config);
                break;
            }

            default:
                res = create(opt, amplicon_storage, pool_id, threshold, config);

        }

        return res;

    }

    AuxiliaryData* AuxiliaryDataFactory::create(const AuxiliaryDataOption opt, const AmpliconStorage& amplicon_storage,
            const numSeqs_t pool_id, const dist_t threshold, const lenSeqs_t num_extra_segments, const SpaceLevenshteinConfiguration& config) {

        typedef BitVector<size_t> bit_vector_t;
        typedef BitRankOne<size_t, 20> rank_t;

        AuxiliaryData* res = nullptr;

        switch (opt) {

            // Reference:
            // Auxiliary data structures employing a (one-way) segment filter to efficiently search
            // for (not yet swarmed) partners of an amplicon.
            case AD_SEGMENT_FILTER: {
                res = new SegmentFilterAuxiliaryData(amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                        num_extra_segments, config);
                break;
            }

            /* The remaining cases use space-efficient variants of the segment filter. */

            // Auxiliary data structures employing a (one-way) segment filter based on NaiveK2Tree instances.
            case AD_SPACE_SEGMENT_FILTER_NAIVE: {
                if (config.tree_parameters.static_rels) {
                    res = new SpaceSegmentFilterAuxiliaryData<StaticSpaceSegmentRelations<NaiveK2Tree<numSeqs_t>>>(
                        amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                        num_extra_segments, config);
                } else {
                    res = new SpaceSegmentFilterAuxiliaryData<SpaceSegmentRelations<NaiveK2Tree<numSeqs_t>>>(
                        amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                        num_extra_segments, config);
                }

                break;
            }

            // Auxiliary data structures employing a (one-way) segment filter based on BasicK2Tree instances.
            case AD_SPACE_SEGMENT_FILTER_BASIC: {
                if (config.tree_parameters.static_rels) {
                    res = new SpaceSegmentFilterAuxiliaryData<StaticSpaceSegmentRelations<BasicK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                            amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                            num_extra_segments, config);
                } else {
                    res = new SpaceSegmentFilterAuxiliaryData<SpaceSegmentRelations<BasicK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                            amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                            num_extra_segments, config);
                }
                break;
            }

            // Auxiliary data structures employing a (one-way) segment filter based on HybridK2Tree instances.
            case AD_SPACE_SEGMENT_FILTER_HYBRID: {
                if (config.tree_parameters.static_rels) {
                    res = new SpaceSegmentFilterAuxiliaryData<StaticSpaceSegmentRelations<HybridK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                            amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                            num_extra_segments, config);
                } else {
                    res = new SpaceSegmentFilterAuxiliaryData<SpaceSegmentRelations<HybridK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                            amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                            num_extra_segments, config);
                }
                break;
            }

            // Auxiliary data structures employing a (one-way) segment filter based on RectangularK2Tree instances.
            case AD_SPACE_SEGMENT_FILTER_RECT: {
                if (config.tree_parameters.static_rels) {
                    res = new SpaceSegmentFilterAuxiliaryData<StaticSpaceSegmentRelations<RectangularK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                            amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                            num_extra_segments, config);
                } else {
                    res = new SpaceSegmentFilterAuxiliaryData<SpaceSegmentRelations<RectangularK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                            amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                            num_extra_segments, config);
                }
                break;
            }

            // Auxiliary data structures employing a (one-way) segment filter based on PartitionedRectangularK2Tree instances.
            case AD_SPACE_SEGMENT_FILTER_PART: {
                if (config.tree_parameters.static_rels) {
                    res = new SpaceSegmentFilterAuxiliaryData<StaticSpaceSegmentRelations<PartitionedRectangularK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                            amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                            num_extra_segments, config);
                } else {
                    res = new SpaceSegmentFilterAuxiliaryData<SpaceSegmentRelations<PartitionedRectangularK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                            amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                            num_extra_segments, config);
                }
                break;
            }

            // Auxiliary data structures employing a (one-way) segment filter based on PartitionedVariableK2Tree instances.
            case AD_SPACE_SEGMENT_FILTER_VARI: {
                if (config.tree_parameters.static_rels) {
                    res = new SpaceSegmentFilterAuxiliaryData<StaticSpaceSegmentRelations<PartitionedVariableK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                            amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                            num_extra_segments, config);
                } else {
                    res = new SpaceSegmentFilterAuxiliaryData<SpaceSegmentRelations<PartitionedVariableK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                            amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                            num_extra_segments, config);
                }
                break;
            }

            default:
                res = create(opt, amplicon_storage, pool_id, threshold, config);

        }

        return res;

    }

    AuxiliaryData* AuxiliaryDataFactory::create(const AuxiliaryDataOption opt, const AmpliconStorage& amplicon_storage,
                                                const numSeqs_t pool_id, const dist_t threshold, const AlignmentFreeConfiguration& config) {

        AuxiliaryData* res = nullptr;

        switch (opt) {

            // Auxiliary data structures providing the bare minimum to execute the clustering phase.
            // No filtering of potential partners.
            case AD_NAIVE_AUXILIARY: {
                res = new FeatureAuxiliaryData(dynamic_cast<const FeatureAmpliconStorage&>(amplicon_storage), pool_id, threshold, config);
                break;
            }

            // Auxiliary data structures employing a space-partitioning structure to efficiently search
            // for (not yet swarmed) partners of an amplicon using the feature representations of the amplicons.
            // Can only be used with norm-like distance functions.
            case AD_KDTREE_AUXILIARY: {
                auto& fas = dynamic_cast<const FeatureAmpliconStorage&>(amplicon_storage);
                auto num_features = fas.get_pool(pool_id).num_features();

                switch (config.get_opt_distance()) {

                    case DT_MANHATTAN: {
                        res = new SpacePartitioningAuxiliaryData<KdTree<manhattan_norm_tag>>(
                                fas, pool_id, threshold, config);
                        break;
                    }

                    case DT_EUCLIDEAN: {
                        res = new SpacePartitioningAuxiliaryData<KdTree<euclidean_norm_tag>>(
                                fas, pool_id, threshold, config);
                        break;
                    }

                    default:
                        // should be caught by check(...) of AlignmentFreeConfiguration
                        std::cerr << "ERROR: Selected 'distance' option is not compatible with the space-partitioning structure." << std::endl;

                }
                break;
            }

            case AD_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown AuxiliaryData option." << std::endl;

        }

        return res;

    }

    AuxiliaryData* AuxiliaryDataFactory::create(const RefinementAuxiliaryDataOption opt, const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, const numSeqs_t pool_id, const dist_t threshold, const Configuration& config) {

        AuxiliaryData* res = nullptr;

        switch (opt) {

            // Auxiliary data structures providing the bare minimum to execute the refinement phase.
            // No filtering of potential partners.
            case RD_NAIVE_AUXILIARY: {
                res = new SimpleAuxiliaryData(amplicon_storage, swarm_storage, pool_id, threshold, config);
                break;
            }

            case RD_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown AuxiliaryData option." << std::endl;

        }

        return res;

    }

    AuxiliaryData* AuxiliaryDataFactory::create(const RefinementAuxiliaryDataOption opt, const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, const numSeqs_t pool_id, const dist_t threshold, const lenSeqs_t num_extra_segments,
            const Configuration& config) {

        AuxiliaryData* res = nullptr;

        switch (opt) {

            // Auxiliary data structures employing a (one-way) segment filter to efficiently search
            // for potential grafting links.
            case RD_SEGMENT_FILTER: {
                res = new SegmentFilterAuxiliaryData(amplicon_storage, swarm_storage, pool_id,
                        static_cast<lenSeqs_t>(threshold), num_extra_segments, config);
                break;
            }

            // Auxiliary data structures employing a (two-way) segment filter to efficiently search
            // for potential grafting links.
            case RD_2W_SEGMENT_FILTER: {
                res = new TwoWaySegmentFilterAuxiliaryData(amplicon_storage, swarm_storage, pool_id,
                        static_cast<lenSeqs_t>(threshold), num_extra_segments, config);
                break;
            }

            // Auxiliary data structures employing a (one-way) segment filter to efficiently search
            // for potential grafting links.
            // The set-up process of the segment filter is adapted in order to determine the number of segments
            // needed for the filter based on the distance threshold and the edit-operation costs.
            case RD_SCORE_SEGMENT_FILTER: {
                ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
                lenSeqs_t num_threshold_segments = static_cast<lenSeqs_t>(threshold) / std::min({sf.penalty_mismatch, sf.penalty_open, sf.penalty_extend});

                res = new ScoreSegmentFilterAuxiliaryData(amplicon_storage, swarm_storage, pool_id,
                        static_cast<lenSeqs_t>(threshold), num_threshold_segments, num_extra_segments, config);
                break;
            }

            default:
                res = create(opt, amplicon_storage, swarm_storage, pool_id, threshold, config);

        }

        return res;

    }

    AuxiliaryData* AuxiliaryDataFactory::create(const RefinementAuxiliaryDataOption opt, const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, const numSeqs_t pool_id, const dist_t threshold, const lenSeqs_t num_extra_segments,
            const SpaceLevenshteinConfiguration& config) {

        typedef BitVector<size_t> bit_vector_t;
        typedef BitRankOne<size_t, 20> rank_t;

        AuxiliaryData* res = nullptr;

        switch (opt) {

            // Reference:
            // Auxiliary data structures employing a (one-way) segment filter to efficiently search
            // for potential grafting links.
            case RD_SEGMENT_FILTER: {
                res = new SegmentFilterAuxiliaryData(amplicon_storage, swarm_storage, pool_id,
                        static_cast<lenSeqs_t>(threshold), num_extra_segments, config);
                break;
            }

            /* The remaining cases use space-efficient variants of the segment filter. */

            // Auxiliary data structures employing a (one-way) segment filter based on NaiveK2Tree instances.
            case RD_SPACE_SEGMENT_FILTER_NAIVE: {
                res = new SpaceSegmentFilterAuxiliaryData<SpaceSegmentRelations<NaiveK2Tree<numSeqs_t>>>(
                        amplicon_storage, swarm_storage, pool_id,
                        static_cast<lenSeqs_t>(threshold), num_extra_segments, config);
                break;
            }

            // Auxiliary data structures employing a (one-way) segment filter based on BasicK2Tree instances.
            case RD_SPACE_SEGMENT_FILTER_BASIC: {
                res = new SpaceSegmentFilterAuxiliaryData<SpaceSegmentRelations<BasicK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                        amplicon_storage, swarm_storage, pool_id,
                        static_cast<lenSeqs_t>(threshold), num_extra_segments, config);
                break;
            }

            // Auxiliary data structures employing a (one-way) segment filter based on HybridK2Tree instances.
            case RD_SPACE_SEGMENT_FILTER_HYBRID: {
                res = new SpaceSegmentFilterAuxiliaryData<SpaceSegmentRelations<HybridK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                        amplicon_storage, swarm_storage, pool_id,
                        static_cast<lenSeqs_t>(threshold), num_extra_segments, config);
                break;
            }

            // Auxiliary data structures employing a (one-way) segment filter based on RectangularK2Tree instances.
            case RD_SPACE_SEGMENT_FILTER_RECT: {
                res = new SpaceSegmentFilterAuxiliaryData<SpaceSegmentRelations<RectangularK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                        amplicon_storage, swarm_storage, pool_id,
                        static_cast<lenSeqs_t>(threshold), num_extra_segments, config);
                break;
            }

            // Auxiliary data structures employing a (one-way) segment filter based on PartitionedRectangularK2Tree instances.
            case RD_SPACE_SEGMENT_FILTER_PART: {
                res = new SpaceSegmentFilterAuxiliaryData<SpaceSegmentRelations<PartitionedRectangularK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                        amplicon_storage, swarm_storage, pool_id,
                        static_cast<lenSeqs_t>(threshold), num_extra_segments, config);
                break;
            }

            // Auxiliary data structures employing a (one-way) segment filter based on PartitionedVariableK2Tree instances.
            case RD_SPACE_SEGMENT_FILTER_VARI: {
                res = new SpaceSegmentFilterAuxiliaryData<SpaceSegmentRelations<PartitionedVariableK2Tree<bit_vector_t, rank_t, numSeqs_t>>>(
                        amplicon_storage, swarm_storage, pool_id,
                        static_cast<lenSeqs_t>(threshold), num_extra_segments, config);
                break;
            }

            default:
                res = create(opt, amplicon_storage, swarm_storage, pool_id, threshold, config);

        }

        return res;

    }

    AuxiliaryData* AuxiliaryDataFactory::create(const RefinementAuxiliaryDataOption opt, const AmpliconStorage& amplicon_storage,
                                                const SwarmStorage& swarm_storage, const numSeqs_t pool_id, const dist_t threshold,
                                                const AlignmentFreeConfiguration& config) {

        AuxiliaryData* res = nullptr;

        switch (opt) {

            // Auxiliary data structures providing the bare minimum to execute the refinement phase.
            // No filtering of potential partners.
            case RD_NAIVE_AUXILIARY: {
                res = new FeatureAuxiliaryData(dynamic_cast<const FeatureAmpliconStorage&>(amplicon_storage), swarm_storage, pool_id, threshold, config);
                break;
            }

            // Auxiliary data structures employing a space-partitioning structure to efficiently search
            // for (not yet swarmed) partners of an amplicon using the feature representations of the amplicons.
            // Can only be used with norm-like distance functions.
            case RD_KDTREE_AUXILIARY: {
                auto& fas = dynamic_cast<const FeatureAmpliconStorage&>(amplicon_storage);
                auto num_features = fas.get_pool(pool_id).num_features();

                switch (config.get_opt_distance()) {

                    case DT_MANHATTAN: {
                        res = new SpacePartitioningAuxiliaryData<KdTree<manhattan_norm_tag>>(
                                fas, swarm_storage, pool_id, threshold, config);
                        break;
                    }

                    case DT_EUCLIDEAN: {
                        res = new SpacePartitioningAuxiliaryData<KdTree<euclidean_norm_tag>>(
                                fas, swarm_storage, pool_id, threshold, config);
                        break;
                    }

                    default:
                        // should be caught by check(...) of AlignmentFreeConfiguration
                        std::cerr << "ERROR: Selected 'distance' option is not compatible with the space-partitioning structure." << std::endl;

                }
                break;
            }

            case RD_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown AuxiliaryData option." << std::endl;

        }

        return res;

    }

}



