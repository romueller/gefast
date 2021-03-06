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

#ifndef GEFAST_QUALITYLEVENSHTEINMODE_HPP
#define GEFAST_QUALITYLEVENSHTEINMODE_HPP

#include "../Base.hpp"

namespace GeFaST {

    /*
     * Configuration of GeFaST's quality-weighted Levenshtein mode (qlev).
     *
     * Computes the distance between amplicons as the number of edit operations in an quality-weighted,
     * optimal alignment using a given scoring function.
     * Uses a segment filter (and, optionally, a q-gram filter) to speed up the clustering phases.
     * Several quality-weighting mechanisms and boosting functions are available.
     */
    struct QualityLevenshteinConfiguration : public Configuration {

        QualityLevenshteinConfiguration(int argc, const char* argv[]);

        QualityLevenshteinConfiguration* clone() const override; // deep-copy clone method


        Preprocessor* build_preprocessor() const override;

        QualityEncoding<>* build_quality_encoding() const override;

        Clusterer* build_clusterer() const override;

        ClusterRefiner* build_cluster_refiner() const override;

        OutputGenerator* build_output_generator() const override;

        AmpliconStorage* build_amplicon_storage(const DataStatistics<>& ds) const override;

        SwarmStorage* build_swarm_storage(const AmpliconStorage& amplicon_storage) const override;

        AuxiliaryData* build_auxiliary_data(const AmpliconStorage& amplicon_storage, const numSeqs_t pool_id,
                const dist_t threshold) const override;

        AuxiliaryData* build_auxiliary_data(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                const numSeqs_t pool_id, const dist_t threshold) const override;

        Distance* build_distance_function(const AmpliconStorage &amplicon_storage, const dist_t threshold) const override;


        void print(std::ostream& stream) const override;

        /*
         * Check whether the configuration obeys the restrictions of the mode.
         */
        void check() const;


        /* Mode-specific configuration parameters */

        lenSeqs_t num_extra_segments; // parameter of the pigeonhole principle (segment filter)
        bool two_way_segment_filter; // mode of the segment filter (false = forward, true = forward-backward)
        bool use_qgrams; // flag indicating whether alignment computations are preceded with a q-gram-based filtering step

        QualityEncoding<> qe_map; // expected quality encoding

        // quality-weighting and boosting
        bool unweighted_matches; // Boolean flag indicating whether match operations are quality-weighted (true iff they remain unweighted)
        bool inner_boost; // Boolean flag indicating whether inner or outer boosting is used (true iff inner boosting is requested)
        BoostingOption opt_boosting; // requested boosting function

    };

}

#endif //GEFAST_QUALITYLEVENSHTEINMODE_HPP
