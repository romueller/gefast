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

#include <iostream>

#include "../include/Factories.hpp"
#include "../include/modes/DereplicationMode.hpp"
#include "../include/AmpliconCollections.hpp"
#include "../include/AmpliconStorages.hpp"
#include "../include/Distances.hpp"
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

            case CF_CONSISTENCY:
                config = new ConsistencyConfiguration(argc, argv);
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



    Preprocessor* PreprocessorFactory::create(const PreprocessorOption opt, const QualityEncoding<>& qe,
            const Configuration& config) {

        Preprocessor* res = nullptr;

        switch (opt) {
            case PP_FASTA:
                res = new FastaPreprocessor();
                break;

            case PP_FASTQ:
                res = new FastqPreprocessor(qe);
                break;

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

            case PP_FASTQ:
                res = new FastqPreprocessor(opt_qe);
                break;

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
                std::cerr << "ERROR: Unknown ClusterRefiner option." << std::endl;

        }

        return res;

    }



    SwarmStorage* SwarmStorageFactory::create(const SwarmStorageOption opt, const AmpliconStorage& amplicon_storage,
            const Configuration& config) {

        SwarmStorage* res = nullptr;

        switch (opt) {

            case SS_SIMPLE_PER_POOL: {
                res = new SimpleSwarmStorage(amplicon_storage.num_pools());
                break;
            }

            case SS_UNKNOWN: // fallthrough to default for error handling
            default:
                std::cerr << "ERROR: Unknown ClusterRefiner option." << std::endl;

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
                std::cerr << "ERROR: Unknown ClusterRefiner option." << std::endl;

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
                res = new ScoreSegmentFilterAuxiliaryData(amplicon_storage, pool_id, static_cast<lenSeqs_t>(threshold),
                        num_extra_segments, config);
                break;
            }

            default:
                res = create(opt, amplicon_storage, pool_id, threshold, config);

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
                std::cerr << "ERROR: Unknown ClusterRefiner option." << std::endl;

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
                res = new ScoreSegmentFilterAuxiliaryData(amplicon_storage, swarm_storage, pool_id,
                        static_cast<lenSeqs_t>(threshold), num_extra_segments, config);
                break;
            }

            default:
                res = create(opt, amplicon_storage, swarm_storage, pool_id, threshold, config);

        }

        return res;

    }

}



