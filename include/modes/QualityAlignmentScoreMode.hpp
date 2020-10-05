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

#ifndef GEFAST_QUALITYALIGNMENTSCOREMODE_HPP
#define GEFAST_QUALITYALIGNMENTSCOREMODE_HPP

#include "../Base.hpp"

namespace GeFaST {

    /*
     * Configuration of GeFaST's quality-weighted alignment-score mode (qas).
     *
     * Computes the distance between amplicons as the (bounded) score of a quality-weighted, optimal alignment between them.
     * The number of bands (diagonals of the dynamic-programming matrix) can be limited to speed up the computations.
     * Also, a segment filter similar to the Levenshtein mode can be activated.
     * However, the segment filter is not affected by the quality-weighting and, thus, rejects amplicon pairs
     * with an unweighted score / distance above the threshold (even though the weighted one is below the threshold).
     * Several quality-weighting mechanisms and boosting functions are available.
     */
    struct QualityAlignmentScoreConfiguration : public Configuration {

        QualityAlignmentScoreConfiguration(int argc, const char* argv[]);

        QualityAlignmentScoreConfiguration* clone() const override; // deep-copy clone method


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

        long long bands_per_side; // number of bands (diagonals) considered on each side of the main diagonal in the DP-matrix
        lenSeqs_t num_extra_segments; // parameter of the pigeonhole principle (segment filter; 0 when the segment filter is not used)
        bool use_qgrams; // flag indicating whether alignment computations are preceded with a q-gram-based filtering step

        QualityEncoding<> qe_map; // expected quality encoding

        // quality-weighting and boosting
        bool unweighted_matches; // Boolean flag indicating whether match operations are quality-weighted (true iff they remain unweighted)
        bool inner_boost; // Boolean flag indicating whether inner or outer boosting is used (true iff inner boosting is requested)
        BoostingOption opt_boosting; // requested boosting function

    };

}

#endif //GEFAST_QUALITYALIGNMENTSCOREMODE_HPP
