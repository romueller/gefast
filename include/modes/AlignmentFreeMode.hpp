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

#ifndef GEFAST_ALIGNMENTFREEMODE_HPP
#define GEFAST_ALIGNMENTFREEMODE_HPP

#include "../Base.hpp"
#include "../Distances.hpp"

namespace GeFaST {

    /*
     * Configuration of GeFaST's alignment-free mode (af).
     *
     * Computes the distance between amplicons as the distance between their features.
     * Several functions for the distance are available.
     */
    struct AlignmentFreeConfiguration : public Configuration {


        AlignmentFreeConfiguration(int argc, const char* argv[]);

        AlignmentFreeConfiguration* clone() const override; // deep-copy clone method


        Preprocessor* build_preprocessor() const override;

        QualityEncoding<>* build_quality_encoding() const override;

        Clusterer* build_clusterer() const override;

        ClusterRefiner* build_cluster_refiner() const override;

        OutputGenerator* build_output_generator() const override;

        // covariant return type compared to Configuration
        FeatureAmpliconStorage* build_amplicon_storage(const DataStatistics<>& ds) const override;

        SwarmStorage* build_swarm_storage(const AmpliconStorage& amplicon_storage) const override;

        AuxiliaryData* build_auxiliary_data(const AmpliconStorage& amplicon_storage, const numSeqs_t pool_id,
                                            const dist_t threshold) const override;

        AuxiliaryData* build_auxiliary_data(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                                            const numSeqs_t pool_id, const dist_t threshold) const override;

        // covariant return type compared to Configuration
        FeatureDistance* build_distance_function(const AmpliconStorage &amplicon_storage, const dist_t threshold) const override;


        void print(std::ostream& stream) const override;

        DistanceOption get_opt_distance() const;


        /* Mode-specific configuration parameters */

        // currently none (except opt_feature, see below)


    protected:
        /*
         * Check whether the configuration obeys the restrictions of the mode.
         */
        void check(int argc, const char* argv[]) const;


        // additional data structure options
        FeatureBuilderOption opt_feature; // features used to represent amplicons

    };

}

#endif //GEFAST_ALIGNMENTFREEMODE_HPP
